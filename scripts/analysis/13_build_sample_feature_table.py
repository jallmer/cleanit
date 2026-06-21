#!/usr/bin/env python3
"""
Build a per-SRR feature table for sample-specific trimming analysis.
"""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path


OUT_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis")
COUNTS_BASE = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/flattened_counts")
TRIM_FILE = OUT_DIR / "trimmomatic_detail.tsv"
BOWTIE_FILE = OUT_DIR / "bowtie2_alignment_stats.tsv"
QUALITY_FILE = OUT_DIR / "per_srr_quality.tsv"
BENEFIT_FILE = OUT_DIR / "trimming_benefit.tsv"
CLASS_FILE = OUT_DIR / "trimming_classification.tsv"
OUT_FILE = OUT_DIR / "sample_feature_table.tsv"

METHOD_MAP = {"U": "untrimmed", "A": "adapter_only", "P5": "P5", "P10": "P10", "P20": "P20", "P35": "P35"}
MODE_PATTERNS = {
    "untrimmed": "untrmd_{srr}_fC.txt.gz",
    "adapter_only": "{srr}_trimmomatic_adapter_fC.txt.gz",
    "P5": "{srr}_trimmomatic_P5_fC.txt.gz",
    "P10": "{srr}_trimmomatic_P10_fC.txt.gz",
    "P20": "{srr}_trimmomatic_P20_fC.txt.gz",
    "P35": "{srr}_trimmomatic_P35_fC.txt.gz",
}


def load_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def count_detected_genes(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    detected = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 7 and int(parts[6]) > 0:
                detected += 1
    return detected


def load_summary_metrics(project: str, srr: str, mode: str) -> tuple[float | None, int | None]:
    base = COUNTS_BASE / project
    stem = MODE_PATTERNS[mode].format(srr=srr)
    summary_paths = [base / stem.replace(".txt.gz", ".txt.summary.gz"), base / stem.replace(".txt.gz", ".txt.summary")]
    assigned = total = None
    for path in summary_paths:
        if path.exists():
            opener = gzip.open if path.suffix == ".gz" else open
            assigned = 0
            total = 0
            with opener(path, "rt") as fh:
                next(fh, None)
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        val = int(parts[1])
                        total += val
                        if parts[0] == "Assigned":
                            assigned = val
            break
    count_path = base / stem
    detected = count_detected_genes(count_path) if count_path.exists() else None
    assigned_frac = (assigned / total) if assigned is not None and total else None
    return assigned_frac, detected


def main() -> int:
    quality_rows = load_tsv(QUALITY_FILE)
    benefit_rows = load_tsv(BENEFIT_FILE)
    class_rows = load_tsv(CLASS_FILE)
    trim_rows = load_tsv(TRIM_FILE)
    bowtie_rows = load_tsv(BOWTIE_FILE)

    quality_map = {r["SRR_ID"]: r for r in quality_rows}
    benefit_map = {r["SRR_ID"]: r for r in benefit_rows}
    class_map = {r["SRR_ID"]: r for r in class_rows}

    trim_map = {}
    for row in trim_rows:
        trim_map[(row["SRR_ID"], row["mode"])] = row

    bowtie_map = {}
    for row in bowtie_rows:
        bowtie_map[(row["srr_id"], row["mode"])] = row

    fieldnames = [
        "SRR_ID", "project_id", "t_star", "benefit_B", "benefit_B_net",
        "Q_mean", "Q_median", "sequence_depth", "read_length_mean",
        "frac_below_q20", "frac_below_q30", "tail_quality_decay",
        "adapter_rate", "duplication_rate", "n_content", "gc_content", "gc_deviation",
    ]
    for method, mode in METHOD_MAP.items():
        fieldnames.extend([
            f"{method}_read_loss",
            f"{method}_assigned_frac",
            f"{method}_detected_genes",
            f"{method}_overall_alignment_rate",
            f"{method}_class",
            f"{method}_delta_path",
        ])

    all_srrs = sorted(set(quality_map) | set(benefit_map) | set(class_map))
    tmp_file = OUT_FILE.with_suffix(".tmp")
    with tmp_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, srr in enumerate(all_srrs, start=1):
            q = quality_map.get(srr, {})
            b = benefit_map.get(srr, {})
            c = class_map.get(srr, {})
            project = b.get("project_id") or q.get("project_id") or c.get("project_id") or ""
            row = {
                "SRR_ID": srr,
                "project_id": project,
                "t_star": c.get("t_star", "NA"),
                "benefit_B": b.get("benefit_B", "NA"),
                "benefit_B_net": b.get("benefit_B_net", "NA"),
            }
            for key in ["Q_mean", "Q_median", "sequence_depth", "read_length_mean", "frac_below_q20", "frac_below_q30",
                        "tail_quality_decay", "adapter_rate", "duplication_rate", "n_content", "gc_content", "gc_deviation"]:
                row[key] = q.get(key, "NA")

            for method, mode in METHOD_MAP.items():
                trim = trim_map.get((srr, mode), {})
                
                bow_key = None
                if method == "U":
                    if (srr, "untrimmed") in bowtie_map:
                        bow_key = "untrimmed"
                    elif (srr, f"untrmd_{srr}") in bowtie_map:
                        bow_key = f"untrmd_{srr}"
                else:
                    mapping = {"A": "adapter"}
                    sfx = mapping.get(method, method)
                    bow_key = f"{srr}_trimmomatic_{sfx}"
                
                bow = bowtie_map.get((srr, bow_key), {})
                assigned_frac, detected = load_summary_metrics(project, srr, mode) if project else (None, None)
                try:
                    input_reads = float(trim.get("input_reads", ""))
                    surviving = float(trim.get("surviving", ""))
                    read_loss = 1.0 - (surviving / input_reads) if input_reads else None
                except Exception:
                    read_loss = None
                row[f"{method}_read_loss"] = f"{read_loss:.6f}" if read_loss is not None else "NA"
                row[f"{method}_assigned_frac"] = f"{assigned_frac:.6f}" if assigned_frac is not None else "NA"
                row[f"{method}_detected_genes"] = str(detected) if detected is not None else "NA"
                rate = bow.get("overall_alignment_rate", "").rstrip("%")
                row[f"{method}_overall_alignment_rate"] = rate if rate else "NA"
                row[f"{method}_class"] = c.get(f"{method}_class", "NA")
                row[f"{method}_delta_path"] = c.get(f"{method}_delta_path", "NA")

            writer.writerow(row)
            if idx % 50 == 0:
                print(f"  Processed {idx}/{len(all_srrs)} SRRs...", flush=True)

    os.replace(tmp_file, OUT_FILE)
    print(f"Written: {OUT_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
