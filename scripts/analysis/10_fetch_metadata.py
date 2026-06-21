#!/usr/bin/env python3
"""
10_fetch_metadata.py — Fetch sample condition metadata from NCBI/ENA for BioProjects.

Queries the NCBI Entrez API to retrieve biosample attributes and attempts
to extract condition/treatment groupings for DE analysis.

Usage:
    python3 scripts/analysis/10_fetch_metadata.py [--project PRJNA...]

Output is written to ~/scripts/analysis/metadata/<PROJECT>.csv with columns:
    Run,condition

These files must be reviewed manually before use in the DE pipeline,
as automated extraction may not correctly identify the relevant contrast.
"""

import os
import sys
import csv
import json
import time
import argparse
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

METADATA_DIR = os.path.join(os.path.expanduser("~"), "scripts", "analysis", "metadata")
PROJECT_IDS_FILE = os.path.expanduser("~/scripts/analysis/project_ids.txt")

# NCBI Entrez base URLs
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ENA API
ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Rate limiting
RATE_LIMIT_SECONDS = 0.4


def fetch_url(url, retries=3):
    """Fetch URL with retries."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cleanit-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"  WARN: Failed to fetch {url}: {e}", file=sys.stderr)
                return None
    return None


def fetch_ena_metadata(project_id):
    """
    Fetch sample metadata from ENA for a BioProject.
    Returns: list of dicts with 'run_accession' and various attributes.
    """
    params = urllib.parse.urlencode({
        "accession": project_id,
        "result": "read_run",
        "fields": "run_accession,sample_accession,sample_title,experiment_title,library_strategy,library_source",
        "format": "tsv",
        "limit": "0",
    })
    url = f"{ENA_FILEREPORT_URL}?{params}"
    content = fetch_url(url)
    if not content:
        return []

    rows = []
    lines = content.strip().split("\n")
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= len(header):
            rows.append(dict(zip(header, parts)))
    return rows


def fetch_ncbi_biosample_attributes(biosample_ids):
    """
    Fetch BioSample attributes from NCBI for a list of biosample IDs.
    Returns: dict biosample_id -> {attribute_name: value}
    """
    results = {}
    batch_size = 50

    for i in range(0, len(biosample_ids), batch_size):
        batch = biosample_ids[i:i + batch_size]
        ids_str = ",".join(batch)
        params = urllib.parse.urlencode({
            "db": "biosample",
            "id": ids_str,
            "rettype": "xml",
            "retmode": "xml",
        })
        url = f"{EFETCH_URL}?{params}"
        content = fetch_url(url)
        if not content:
            continue

        try:
            root = ET.fromstring(content)
            for sample in root.findall(".//BioSample"):
                accession = sample.get("accession", "")
                attrs = {}
                for attr in sample.findall(".//Attribute"):
                    name = attr.get("attribute_name", attr.get("harmonized_name", ""))
                    value = attr.text or ""
                    if name:
                        attrs[name.lower().replace(" ", "_")] = value.strip()
                # Also get the title
                title_elem = sample.find(".//Title")
                if title_elem is not None and title_elem.text:
                    attrs["title"] = title_elem.text.strip()
                results[accession] = attrs
        except ET.ParseError as e:
            print(f"  WARN: XML parse error: {e}", file=sys.stderr)

        time.sleep(RATE_LIMIT_SECONDS)

    return results


def infer_condition(attrs):
    """
    Attempt to infer a condition/group label from BioSample attributes.
    Returns the best candidate condition string.
    """
    # Priority order for condition attributes
    condition_keys = [
        "treatment", "condition", "disease", "disease_state",
        "phenotype", "genotype", "cell_type", "tissue",
        "source_name", "sample_type", "group", "age",
        "histology", "diagnosis", "clinical_stage",
    ]

    for key in condition_keys:
        if key in attrs and attrs[key]:
            val = attrs[key].strip()
            if val.lower() not in ("not applicable", "n/a", "na", "missing", "unknown", ""):
                return val

    # Fallback: use title
    if "title" in attrs:
        return attrs["title"]

    return "unknown"


def process_project(project_id):
    """
    Fetch metadata for a BioProject and write condition CSV.
    """
    print(f"  Fetching ENA metadata for {project_id}...")
    ena_rows = fetch_ena_metadata(project_id)
    if not ena_rows:
        print(f"  No ENA data found for {project_id}")
        return False

    print(f"  Found {len(ena_rows)} runs")

    # Extract unique biosample accessions
    biosample_ids = list(set(r.get("sample_accession", "") for r in ena_rows if r.get("sample_accession")))

    if biosample_ids:
        print(f"  Fetching BioSample attributes for {len(biosample_ids)} samples...")
        biosample_attrs = fetch_ncbi_biosample_attributes(biosample_ids)
    else:
        biosample_attrs = {}

    # Map run -> condition
    run_conditions = []
    for row in ena_rows:
        run = row.get("run_accession", "")
        sample = row.get("sample_accession", "")
        attrs = biosample_attrs.get(sample, {})

        # Add sample title from ENA as fallback
        if "sample_title" in row and row["sample_title"]:
            attrs.setdefault("title", row["sample_title"])

        condition = infer_condition(attrs)
        run_conditions.append((run, condition))

    if not run_conditions:
        return False

    # Check if we got meaningful conditions (not all the same)
    conditions = set(c for _, c in run_conditions)
    if len(conditions) <= 1:
        print(f"  WARNING: All samples have the same condition: {conditions}")
        print(f"  This project likely needs manual curation.")

    # Write output
    out_file = os.path.join(METADATA_DIR, f"{project_id}.csv")
    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Run", "condition"])
        for run, cond in sorted(run_conditions):
            writer.writerow([run, cond])
    print(f"  Written: {out_file} ({len(run_conditions)} runs, {len(conditions)} conditions)")

    # Also write full attributes for reference
    attrs_file = os.path.join(METADATA_DIR, f"{project_id}_attributes.tsv")
    if biosample_attrs:
        all_attr_keys = sorted(set(k for attrs in biosample_attrs.values() for k in attrs.keys()))
        with open(attrs_file, "w") as f:
            f.write("biosample_id\t" + "\t".join(all_attr_keys) + "\n")
            for bs_id, attrs in sorted(biosample_attrs.items()):
                f.write(bs_id + "\t" + "\t".join(attrs.get(k, "") for k in all_attr_keys) + "\n")
        print(f"  Full attributes: {attrs_file}")

    return True


def main():
    parser = argparse.ArgumentParser(description="Fetch sample condition metadata from NCBI/ENA")
    parser.add_argument("--project", type=str, default="", help="Fetch only this project")
    parser.add_argument("--project-list", type=str, default="",
                        help="File with one project ID per line")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip projects that already have metadata CSVs")
    args = parser.parse_args()

    os.makedirs(METADATA_DIR, exist_ok=True)

    # Determine which projects to fetch
    projects = []
    if args.project:
        projects = [args.project]
    elif args.project_list:
        with open(args.project_list) as f:
            projects = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    elif os.path.exists(PROJECT_IDS_FILE):
        with open(PROJECT_IDS_FILE) as f:
            projects = [line.strip().split("\t")[0] for line in f if line.strip()]
    else:
        print("No projects specified. Use --project, --project-list, or create project_ids.txt")
        return

    # Skip projects with existing metadata
    if args.skip_existing:
        existing = set()
        for f in os.listdir(METADATA_DIR):
            if f.endswith(".csv") and not f.endswith("_attributes.tsv"):
                existing.add(f.replace(".csv", ""))
        before = len(projects)
        projects = [p for p in projects if p not in existing]
        print(f"Skipping {before - len(projects)} projects with existing metadata")

    print(f"Fetching metadata for {len(projects)} projects")
    print(f"Output directory: {METADATA_DIR}")
    print()

    success = 0
    failed = 0
    for project in projects:
        print(f"[{projects.index(project) + 1}/{len(projects)}] {project}")
        if process_project(project):
            success += 1
        else:
            failed += 1
        time.sleep(RATE_LIMIT_SECONDS)

    print(f"\nDone: {success} succeeded, {failed} failed")
    print(f"\nIMPORTANT: Review the generated CSVs manually before using them.")
    print(f"    The 'condition' column may need curation for correct biological contrast.")


if __name__ == "__main__":
    main()
