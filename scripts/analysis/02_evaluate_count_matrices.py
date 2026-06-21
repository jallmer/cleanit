#!/usr/bin/env python3
"""
02_evaluate_count_matrices.py — Compute JSD and Pearson correlation between count matrices.

For each SRR that has all 6 count files (untrimmed, adapter, P5, P10, P20, P35),
computes pairwise JSD and Pearson correlation.

Usage:
    python3 scripts/analysis/02_evaluate_count_matrices.py [--cores 32] [--project PRJNA...]

Reads from: /scratch/hpc-prf-omiks/ja/flattened_counts/
Writes to:  /scratch/hpc-prf-omiks/ja/analysis/per_srr_eval.tsv
            /scratch/hpc-prf-omiks/ja/analysis/per_project/<PROJECT>_eval.tsv
"""

import os
import sys
import gzip
import glob
import argparse
from collections import defaultdict
from multiprocessing import Pool
import numpy as np

COUNTS_BASE = "/scratch/hpc-prf-omiks/ja/flattened_counts"
OUT_DIR = "/scratch/hpc-prf-omiks/ja/analysis"

MODES = ["untrimmed", "adapter_only", "P5", "P10", "P20", "P35"]

# Map mode names to filename patterns
MODE_PATTERNS = {
    "untrimmed":    "untrmd_{srr}_fC.txt.gz",
    "adapter_only": "{srr}_trimmomatic_adapter_fC.txt.gz",
    "P5":           "{srr}_trimmomatic_P5_fC.txt.gz",
    "P10":          "{srr}_trimmomatic_P10_fC.txt.gz",
    "P20":          "{srr}_trimmomatic_P20_fC.txt.gz",
    "P35":          "{srr}_trimmomatic_P35_fC.txt.gz",
}

HEADER = "\t".join([
    "SRR_ID", "project_id",
    "untrmd_adptrTrmd_jsd", "untrmd_P5Trmd_jsd", "untrmd_P10Trmd_jsd",
    "untrmd_P20Trmd_jsd", "untrmd_P35Trmd_jsd", "P5Trmd_P35Trmd_jsd",
    "untrmd_adptrTrmd_pear", "untrmd_P5Trmd_pear", "untrmd_P10Trmd_pear",
    "untrmd_P20Trmd_pear", "untrmd_P35Trmd_pear", "P5Trmd_P35Trmd_pear",
])


def load_counts(filepath):
    """Load featureCounts output, return dict of gene -> count."""
    counts = {}
    try:
        opener = gzip.open if filepath.endswith(".gz") else open
        with opener(filepath, "rt") as f:
            for line in f:
                if line.startswith("#") or line.startswith("Geneid"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 7:
                    gene = parts[0]
                    count = int(parts[6])  # featureCounts: column 7
                    counts[gene] = count
    except Exception as e:
        print(f"  WARN: Failed to load {filepath}: {e}", file=sys.stderr)
        return None
    return counts


def counts_to_array(counts_a, counts_b):
    """Align two count dicts to matched numpy arrays."""
    genes = sorted(set(counts_a.keys()) | set(counts_b.keys()))
    a = np.array([counts_a.get(g, 0) for g in genes], dtype=np.float64)
    b = np.array([counts_b.get(g, 0) for g in genes], dtype=np.float64)
    return a, b


def jsd(p_counts, q_counts):
    """Jensen-Shannon Divergence between two count vectors (as probability distributions)."""
    p, q = counts_to_array(p_counts, q_counts)
    # Normalize to probability distributions
    p_sum = p.sum()
    q_sum = q.sum()
    if p_sum == 0 or q_sum == 0:
        return np.nan
    p = p / p_sum
    q = q / q_sum
    m = 0.5 * (p + q)
    # KL divergence with epsilon for numerical stability
    eps = 1e-30
    kl_pm = np.sum(p * np.log2((p + eps) / (m + eps)))
    kl_qm = np.sum(q * np.log2((q + eps) / (m + eps)))
    return 0.5 * kl_pm + 0.5 * kl_qm


def pearson_corr(p_counts, q_counts):
    """Pearson correlation between two count vectors."""
    p, q = counts_to_array(p_counts, q_counts)
    if p.std() == 0 or q.std() == 0:
        return np.nan
    return float(np.corrcoef(p, q)[0, 1])


def process_srr(args):
    """Process a single SRR: load 6 count files, compute all metrics."""
    srr, project, counts_dir = args
    
    # Load all 6 count matrices
    data = {}
    for mode in MODES:
        pattern = MODE_PATTERNS[mode].format(srr=srr)
        filepath = os.path.join(counts_dir, pattern)
        if not os.path.exists(filepath):
            return None  # Skip incomplete SRRs
        counts = load_counts(filepath)
        if counts is None:
            return None
        data[mode] = counts

    # Check if P35 has any counts (it might be empty)
    p35_total = sum(data["P35"].values())

    # Compute pairwise comparisons
    # untrimmed vs each mode
    pairs_jsd = {}
    pairs_pear = {}
    for mode in ["adapter_only", "P5", "P10", "P20", "P35"]:
        j = jsd(data["untrimmed"], data[mode])
        p = pearson_corr(data["untrimmed"], data[mode])
        key = f"untrmd_{mode}"
        pairs_jsd[key] = j
        pairs_pear[key] = p

    # P5 vs P35
    j_p5p35 = jsd(data["P5"], data["P35"])
    p_p5p35 = pearson_corr(data["P5"], data["P35"])

    # Format output, using NA for NaN
    def fmt(v):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "NA"
        return str(v)

    row = "\t".join([
        srr, project,
        fmt(pairs_jsd.get("untrmd_adapter_only")),
        fmt(pairs_jsd.get("untrmd_P5")),
        fmt(pairs_jsd.get("untrmd_P10")),
        fmt(pairs_jsd.get("untrmd_P20")),
        fmt(pairs_jsd.get("untrmd_P35")),
        fmt(j_p5p35),
        fmt(pairs_pear.get("untrmd_adapter_only")),
        fmt(pairs_pear.get("untrmd_P5")),
        fmt(pairs_pear.get("untrmd_P10")),
        fmt(pairs_pear.get("untrmd_P20")),
        fmt(pairs_pear.get("untrmd_P35")),
        fmt(p_p5p35),
    ])
    return row


def main():
    parser = argparse.ArgumentParser(description="Evaluate count matrices (JSD + Pearson)")
    parser.add_argument("--cores", type=int, default=32, help="Parallel cores")
    parser.add_argument("--project", type=str, default="", help="Process only this project")
    args = parser.parse_args()

    os.makedirs(os.path.join(OUT_DIR, "per_project"), exist_ok=True)

    # Discover all SRRs across all projects
    work_items = []
    project_dirs = sorted(glob.glob(os.path.join(COUNTS_BASE, "*/")))

    for pdir in project_dirs:
        project = os.path.basename(pdir.rstrip("/"))
        if args.project and project != args.project:
            continue

        # Find unique SRR IDs by looking for untrimmed files
        untrmd_files = glob.glob(os.path.join(pdir, "untrmd_*_fC.txt.gz"))
        for uf in untrmd_files:
            basename = os.path.basename(uf)
            # Pattern: untrmd_SRR12345_fC.txt.gz
            srr = basename.replace("untrmd_", "").replace("_fC.txt.gz", "")
            work_items.append((srr, project, pdir))

    print(f"Found {len(work_items)} SRRs across {len(set(p for _, p, _ in work_items))} projects")
    print(f"Using {args.cores} cores")

    # Process in parallel
    with Pool(processes=args.cores) as pool:
        results = pool.map(process_srr, work_items)

    # Filter out None results (incomplete SRRs)
    valid = [r for r in results if r is not None]
    skipped = len(results) - len(valid)
    print(f"Completed: {len(valid)}, Skipped (incomplete): {skipped}")

    # Write global output
    out_file = os.path.join(OUT_DIR, "per_srr_eval.tsv")
    with open(out_file, "w") as f:
        f.write(HEADER + "\n")
        for row in valid:
            f.write(row + "\n")
    print(f"Written: {out_file}")

    # Write per-project outputs
    project_rows = defaultdict(list)
    for row in valid:
        parts = row.split("\t")
        project_rows[parts[1]].append(row)

    for project, rows in sorted(project_rows.items()):
        pfile = os.path.join(OUT_DIR, "per_project", f"{project}_eval.tsv")
        with open(pfile, "w") as f:
            f.write(HEADER + "\n")
            for row in rows:
                f.write(row + "\n")
    print(f"Written per-project files for {len(project_rows)} projects")


if __name__ == "__main__":
    main()
