#!/usr/bin/env python3
"""
run_bio_concordance.py — Biological concordance analysis across trimming modes.

For each project with ≥2 condition classes and ≥3 samples per class:
  1. Build a count matrix using ALL samples under each trimming mode
  2. Run DESeq2 to get gene-level DE statistics
  3. Run GSEA (prerank) on the Wald statistics
  4. Compare trimmed results against untrimmed reference:
     - Spearman correlation of gene-level statistics
     - Jaccard overlap of significant DEGs
     - Spearman correlation of pathway NES
     - Jaccard overlap of significant pathways
     - Direction concordance of top pathways

Usage:
    python3 run_bio_concordance.py [--project PRJNA...]

Reads from:
    /pc2/users/o/omiks001/hpc-prf-omiks/ja/flattened_counts/
    /pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/deseq2_metadata/sample_sheets/
Writes to:
    /pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/
"""

import os
import sys
import gzip
import csv
import argparse
import warnings
import traceback
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

BASE = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja")
COUNTS_BASE = BASE / "flattened_counts"
METADATA_DIR = BASE / "analysis" / "deseq2_metadata" / "sample_sheets"
OUT_DIR = BASE / "analysis" / "concordance"

MODE_PATTERNS = {
    "U":   "untrmd_{srr}_fC.txt.gz",
    "A":   "{srr}_trimmomatic_adapter_fC.txt.gz",
    "P5":  "{srr}_trimmomatic_P5_fC.txt.gz",
    "P10": "{srr}_trimmomatic_P10_fC.txt.gz",
    "P20": "{srr}_trimmomatic_P20_fC.txt.gz",
    "P35": "{srr}_trimmomatic_P35_fC.txt.gz",
}

ALL_METHODS = ["U", "A", "P5", "P10", "P20", "P35"]
DEG_FDR = 0.05
DEG_LFC = 1.0
GSEA_FDR = 0.25

# Projects that should be excluded or are single-class
EXCLUDE_PROJECTS = {"PRJNA480287", "PRJNA1014106", "PRJNA316201"}

# Only use projects with 2-4 classes for clean biological analysis
MAX_CLASSES = 4
MIN_SAMPLES_PER_CLASS = 3


def load_featurecounts(filepath):
    """Load featureCounts output, return dict gene -> count."""
    counts = {}
    try:
        opener = gzip.open if str(filepath).endswith(".gz") else open
        with opener(filepath, "rt") as f:
            for line in f:
                if line.startswith("#") or line.startswith("Geneid"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 7:
                    counts[parts[0]] = int(parts[6])
    except Exception as e:
        print(f"  WARN: Failed to load {filepath}: {e}", file=sys.stderr)
        return None
    return counts


def load_metadata(project):
    """Load Run,condition CSV for a project."""
    filepath = METADATA_DIR / f"{project}.csv"
    if not filepath.exists():
        return None
    condition_map = {}
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            srr = row["Run"].strip()
            cond = row["condition"].strip()
            condition_map[srr] = cond
    return condition_map


def build_count_matrix(srr_list, project_dir, method, all_genes):
    """
    Build a count matrix DataFrame (genes x samples) for given SRRs and method.
    Returns: (DataFrame[genes x samples], valid_srrs)
    """
    sample_counts = {}
    for srr in srr_list:
        pattern = MODE_PATTERNS[method].format(srr=srr)
        filepath = project_dir / pattern
        if not filepath.exists():
            continue
        counts = load_featurecounts(filepath)
        if counts is None:
            continue
        sample_counts[srr] = counts

    if not sample_counts:
        return None, []

    valid_srrs = list(sample_counts.keys())
    matrix = pd.DataFrame(
        {srr: {g: sample_counts[srr].get(g, 0) for g in all_genes}
         for srr in valid_srrs},
        dtype=np.int64,
    )
    return matrix, valid_srrs


def run_deseq2(count_matrix, conditions, min_total_count=10):
    """
    Run PyDESeq2 on a count matrix.
    count_matrix: DataFrame, genes x samples
    conditions: dict srr -> condition
    Returns: DataFrame with columns [log2FoldChange, stat, pvalue, padj]
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    sample_names = count_matrix.columns.tolist()
    conds = [conditions[s] for s in sample_names]
    unique_conds = sorted(set(conds))
    if len(unique_conds) < 2:
        return None

    # Transpose for PyDESeq2: samples x genes
    count_df = count_matrix.T.copy()
    count_df.index = sample_names

    # Filter low-count genes
    keep = count_df.sum(axis=0) >= min_total_count
    count_df = count_df.loc[:, keep]

    meta_df = pd.DataFrame({"condition": conds}, index=sample_names)
    ref_level = unique_conds[0]

    try:
        dds = DeseqDataSet(
            counts=count_df,
            metadata=meta_df,
            design_factors="condition",
            refit_cooks=True,
            n_cpus=1,
        )
        dds.deseq2()

        stat = DeseqStats(
            dds,
            contrast=["condition", unique_conds[1], ref_level],
            n_cpus=1,
        )
        stat.summary()
        return stat.results_df
    except Exception as e:
        print(f"    DESeq2 failed: {e}", file=sys.stderr)
        return None


HALLMARK_GMT = str(BASE / "analysis" / "MSigDB_Hallmark_2020.gmt")


def run_gsea_prerank(de_results, gene_set_db=None):
    """
    Run GSEA prerank on DE results.
    Returns: DataFrame with columns [Term, NES, FDR q-val, ...]
    """
    import gseapy

    if de_results is None or de_results.empty:
        return None

    if gene_set_db is None:
        gene_set_db = HALLMARK_GMT

    # Use the Wald statistic as ranking metric
    stats = de_results["stat"].dropna()
    stats = stats[stats != 0]
    if len(stats) < 50:
        return None

    rnk = stats.sort_values(ascending=False)

    try:
        pre_res = gseapy.prerank(
            rnk=rnk,
            gene_sets=gene_set_db,
            outdir=None,
            min_size=15,
            max_size=500,
            permutation_num=100,
            seed=42,
            verbose=False,
            no_plot=True,
        )
        return pre_res.res2d
    except Exception as e:
        print(f"    GSEA failed: {e}", file=sys.stderr)
        return None


def jaccard(set_a, set_b):
    """Jaccard index between two sets."""
    if not set_a and not set_b:
        return 1.0
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)


def compute_concordance(ref_de, cand_de, ref_gsea, cand_gsea):
    """
    Compute concordance metrics between reference (untrimmed) and candidate.
    """
    metrics = {}

    # --- Gene-level ---
    common_genes = sorted(set(ref_de.index) & set(cand_de.index))
    if len(common_genes) >= 10:
        ref_stats = ref_de.loc[common_genes, "stat"].fillna(0).values
        cand_stats = cand_de.loc[common_genes, "stat"].fillna(0).values
        rho, _ = spearmanr(ref_stats, cand_stats)
        metrics["rho_gene"] = rho
    else:
        metrics["rho_gene"] = np.nan

    # Significant DEGs
    ref_degs = set(
        ref_de.index[
            (ref_de["padj"] < DEG_FDR) & (ref_de["log2FoldChange"].abs() > DEG_LFC)
        ]
    )
    cand_degs = set(
        cand_de.index[
            (cand_de["padj"] < DEG_FDR) & (cand_de["log2FoldChange"].abs() > DEG_LFC)
        ]
    )
    metrics["n_ref_degs"] = len(ref_degs)
    metrics["n_cand_degs"] = len(cand_degs)
    metrics["jaccard_deg"] = jaccard(ref_degs, cand_degs)

    # --- Pathway-level ---
    if ref_gsea is not None and cand_gsea is not None:
        ref_nes = ref_gsea.set_index("Term")["NES"].astype(float)
        cand_nes = cand_gsea.set_index("Term")["NES"].astype(float)
        common_paths = sorted(set(ref_nes.index) & set(cand_nes.index))

        if len(common_paths) >= 5:
            rho_p, _ = spearmanr(
                ref_nes.loc[common_paths].values,
                cand_nes.loc[common_paths].values,
            )
            metrics["rho_pathway"] = rho_p
        else:
            metrics["rho_pathway"] = np.nan

        # Significant pathways (|NES| > 1 as proxy)
        ref_sig = set(ref_nes.index[ref_nes.abs() > 1.0])
        cand_sig = set(cand_nes.index[cand_nes.abs() > 1.0])
        metrics["jaccard_pathway"] = jaccard(ref_sig, cand_sig)

        # Direction concordance among top 20 pathways
        top_ref = ref_nes.abs().nlargest(20)
        agree = 0
        total = 0
        for path in top_ref.index:
            if path in cand_nes.index:
                total += 1
                if np.sign(ref_nes[path]) == np.sign(cand_nes[path]):
                    agree += 1
        metrics["dir_concordance"] = agree / total if total > 0 else np.nan
        metrics["n_common_pathways"] = len(common_paths)
    else:
        metrics["rho_pathway"] = np.nan
        metrics["jaccard_pathway"] = np.nan
        metrics["dir_concordance"] = np.nan
        metrics["n_common_pathways"] = 0

    return metrics


def process_project(project, condition_map):
    """Process a single project: run DE+GSEA for each method, compute concordance."""
    project_dir = COUNTS_BASE / project
    if not project_dir.exists():
        print(f"  Skipping {project}: count directory not found")
        return None

    # Find SRRs available in both metadata and on disk (check untrimmed as proxy)
    available_srrs = []
    for srr in condition_map:
        fp = project_dir / MODE_PATTERNS["U"].format(srr=srr)
        if fp.exists():
            available_srrs.append(srr)

    if len(available_srrs) < 4:
        print(f"  Skipping {project}: only {len(available_srrs)} samples available")
        return None

    # Check class balance
    cond_counts = pd.Series([condition_map[s] for s in available_srrs]).value_counts()
    if len(cond_counts) < 2:
        print(f"  Skipping {project}: only 1 condition class")
        return None
    if cond_counts.min() < MIN_SAMPLES_PER_CLASS:
        print(f"  Skipping {project}: min class size {cond_counts.min()} < {MIN_SAMPLES_PER_CLASS}")
        return None
    if len(cond_counts) > MAX_CLASSES:
        print(f"  Skipping {project}: {len(cond_counts)} classes > {MAX_CLASSES}")
        return None

    # Only 2-class for this initial pass (cleanest biological comparison)
    conditions = {s: condition_map[s] for s in available_srrs}
    unique_conds = sorted(set(conditions.values()))
    print(f"  {project}: {len(available_srrs)} samples, {len(unique_conds)} classes: {dict(cond_counts)}")

    # Discover gene universe from untrimmed files
    gene_sets = []
    for srr in available_srrs[:3]:  # Sample a few for speed
        fp = project_dir / MODE_PATTERNS["U"].format(srr=srr)
        c = load_featurecounts(fp)
        if c:
            gene_sets.append(set(c.keys()))
    if not gene_sets:
        return None
    all_genes = sorted(set().union(*gene_sets))

    # Run DE+GSEA for each method
    de_results = {}
    gsea_results = {}
    for method in ALL_METHODS:
        print(f"    Method {method}...", end="", flush=True)
        matrix, valid_srrs = build_count_matrix(available_srrs, project_dir, method, all_genes)
        if matrix is None or len(valid_srrs) < 4:
            print(f" skipped (insufficient data)")
            continue

        valid_conditions = {s: conditions[s] for s in valid_srrs}
        de = run_deseq2(matrix, valid_conditions)
        if de is None:
            print(f" DESeq2 failed")
            continue

        gsea = run_gsea_prerank(de)
        de_results[method] = de
        gsea_results[method] = gsea
        n_degs = ((de["padj"] < DEG_FDR) & (de["log2FoldChange"].abs() > DEG_LFC)).sum()
        n_paths = len(gsea) if gsea is not None else 0
        print(f" {n_degs} DEGs, {n_paths} pathways")

    if "U" not in de_results:
        print(f"  Skipping {project}: untrimmed DE failed")
        return None

    # Compute concordance vs untrimmed
    rows = []
    for method in ALL_METHODS:
        if method == "U" or method not in de_results:
            continue
        conc = compute_concordance(
            de_results["U"], de_results[method],
            gsea_results.get("U"), gsea_results.get(method),
        )
        conc["project"] = project
        conc["method"] = method
        conc["n_samples"] = len(available_srrs)
        conc["n_classes"] = len(unique_conds)
        rows.append(conc)

    # Save per-project DE and GSEA results for the notebook
    proj_out = OUT_DIR / project
    proj_out.mkdir(parents=True, exist_ok=True)

    for method, de in de_results.items():
        de.to_csv(proj_out / f"de_{method}.tsv", sep="\t")

    for method, gsea in gsea_results.items():
        if gsea is not None:
            gsea.to_csv(proj_out / f"gsea_{method}.tsv", sep="\t")

    return rows


def main():
    parser = argparse.ArgumentParser(description="Biological concordance analysis")
    parser.add_argument("--project", type=str, default="", help="Process only this project")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Discover eligible projects
    all_results = []
    sheet_files = sorted(METADATA_DIR.glob("*.csv"))
    print(f"Found {len(sheet_files)} sample sheet files")

    for sheet_file in sheet_files:
        project = sheet_file.stem
        if args.project and project != args.project:
            continue
        if project in EXCLUDE_PROJECTS:
            print(f"\n  Skipping {project}: excluded")
            continue

        condition_map = load_metadata(project)
        if condition_map is None:
            continue

        print(f"\n{'='*60}")
        print(f"Project: {project}")
        print(f"{'='*60}")

        try:
            rows = process_project(project, condition_map)
            if rows:
                all_results.extend(rows)
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)

    # Write combined concordance table
    if all_results:
        df = pd.DataFrame(all_results)
        col_order = [
            "project", "method", "n_samples", "n_classes",
            "rho_gene", "jaccard_deg", "n_ref_degs", "n_cand_degs",
            "rho_pathway", "jaccard_pathway", "dir_concordance", "n_common_pathways",
        ]
        col_order = [c for c in col_order if c in df.columns]
        df = df[col_order]
        out_file = OUT_DIR / "bio_concordance.tsv"
        df.to_csv(out_file, sep="\t", index=False)
        print(f"\n{'='*60}")
        print(f"Combined concordance: {out_file} ({len(df)} rows)")
        print(f"{'='*60}")
        print(df.to_string(index=False))
    else:
        print("\nNo results produced.")


if __name__ == "__main__":
    main()
