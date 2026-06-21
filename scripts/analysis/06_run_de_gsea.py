#!/usr/bin/env python3
"""
06_run_de_gsea.py — Leave-one-sample-out DE + GSEA concordance evaluation.

Implements Sections 4–6 of the sample-specific trimming methodology:
  1. For each project with condition metadata, iterate over each target SRR.
  2. Remove the target SRR and run DE + GSEA under each stable reference method.
  3. Compute the consensus reference (median z-scores and NES across stable methods).
  4. Re-insert the target SRR under each candidate method and re-run DE + GSEA.
  5. Compute concordance metrics: Spearman(gene), Spearman(pathway), Jaccard DEG,
     Jaccard pathway, direction concordance.

Usage:
    python3 scripts/analysis/06_run_de_gsea.py [--cores 32] [--project PRJNA...]

Reads from:
    /scratch/hpc-prf-omiks/ja/flattened_counts/
    ~/scripts/analysis/metadata/  (Run,condition CSVs)
Writes to:
    /scratch/hpc-prf-omiks/ja/analysis/concordance/
"""

import os
import sys
import gzip
import glob
import csv
import argparse
import warnings
import traceback
from collections import defaultdict
from multiprocessing import Pool
import numpy as np

# Suppress convergence warnings from PyDESeq2
warnings.filterwarnings("ignore")

COUNTS_BASE = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/flattened_counts"
METADATA_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/deseq2_metadata/sample_sheets"
OUT_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis"
CONCORDANCE_DIR = os.path.join(OUT_DIR, "concordance")
PER_PROJECT_EVAL_DIR = os.path.join(OUT_DIR, "per_project")
LOSO_RAW_DIR = os.path.join(OUT_DIR, "loso_raw_data")
os.makedirs(LOSO_RAW_DIR, exist_ok=True)

# Trimming modes and their filename patterns
MODE_PATTERNS = {
    "U":   "untrmd_{srr}_fC.txt.gz",
    "A":   "{srr}_trimmomatic_adapter_fC.txt.gz",
    "P5":  "{srr}_trimmomatic_P5_fC.txt.gz",
    "P10": "{srr}_trimmomatic_P10_fC.txt.gz",
    "P20": "{srr}_trimmomatic_P20_fC.txt.gz",
    "P35": "{srr}_trimmomatic_P35_fC.txt.gz",
}

# The stable reference methods (Section 4.2)
REF_METHODS = ["U", "A", "P5", "P10"]

# All candidate methods to evaluate (Section 5)
ALL_METHODS = ["U", "A", "P5", "P10", "P20", "P35"]

# Gene set collection for GSEA (Hallmark recommended for main paper)
GENE_SET_DB = "MSigDB_Hallmark_2020"

# Thresholds
DEG_FDR = 0.05
DEG_LFC = 1.0
GSEA_FDR = 0.25
TOP_K_DEG = 200   # Top K genes for Jaccard overlap (Section 6.1b)
TOP_L_PATH = 20   # Top L pathways for direction concordance (Section 6.2e)
DELTA_TOL = 0.01  # Tolerance δ for tie-breaking (Section 8.2)
EPSILON = 0.02    # Threshold ε for helpful/neutral/harmful (Section 9)
MIN_SAMPLES_PER_GROUP = 3
BASELINE_METHOD = "auto"


# ============================================================
# Count Matrix I/O
# ============================================================

def load_featurecounts(filepath):
    """Load featureCounts output, return dict gene -> count."""
    counts = {}
    try:
        opener = gzip.open if filepath.endswith(".gz") else open
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


def build_count_matrix(srr_list, project_dir, method, all_genes=None):
    """
    Build a count matrix DataFrame (genes × samples) for a list of SRRs
    all processed with the same trimming method.
    Returns: (gene_names, count_array[genes × samples], valid_srrs)
    """
    sample_counts = []
    valid_srrs = []

    for srr in srr_list:
        pattern = MODE_PATTERNS[method].format(srr=srr)
        filepath = os.path.join(project_dir, pattern)
        if not os.path.exists(filepath):
            continue
        counts = load_featurecounts(filepath)
        if counts is None:
            continue
        sample_counts.append(counts)
        valid_srrs.append(srr)

    if not sample_counts:
        return None, None, []

    # Union of all genes
    if all_genes is None:
        all_genes = sorted(set().union(*[set(c.keys()) for c in sample_counts]))

    matrix = np.zeros((len(all_genes), len(valid_srrs)), dtype=np.int64)
    for j, counts in enumerate(sample_counts):
        for i, gene in enumerate(all_genes):
            matrix[i, j] = counts.get(gene, 0)

    return all_genes, matrix, valid_srrs


def choose_project_baseline(project_dir, available_srrs, ref_methods):
    """
    Pick the least aggressive method with near-maximal assigned-read fraction.
    This is a lightweight approximation of the project-baseline rule from the methodology.
    """
    method_scores = {}
    order = {"U": 0, "A": 1, "P5": 2, "P10": 3, "P20": 4, "P35": 5}
    for method in ref_methods:
        fractions = []
        for srr in available_srrs:
            stem = MODE_PATTERNS[method].format(srr=srr)
            summary_paths = [
                os.path.join(project_dir, stem.replace(".txt.gz", ".txt.summary.gz")),
                os.path.join(project_dir, stem.replace(".txt.gz", ".txt.summary")),
            ]
            for path in summary_paths:
                if not os.path.exists(path):
                    continue
                opener = gzip.open if path.endswith(".gz") else open
                assigned = 0
                total = 0
                with opener(path, "rt") as f:
                    next(f, None)
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            val = int(parts[1])
                            total += val
                            if parts[0] == "Assigned":
                                assigned = val
                if total > 0:
                    fractions.append(assigned / total)
                break
        if fractions:
            method_scores[method] = float(np.median(fractions))

    if not method_scores:
        return "U"
    best = max(method_scores.values())
    candidates = [m for m, score in method_scores.items() if score >= best - 0.01]
    return min(candidates, key=lambda m: order.get(m, 99))


# ============================================================
# DE Analysis (PyDESeq2)
# ============================================================

def run_deseq2(gene_names, count_matrix, sample_names, conditions):
    """
    Run PyDESeq2 on a count matrix.
    Returns: dict gene -> (log2fc, stat, pvalue, padj)
    """
    try:
        import pandas as pd
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("ERROR: pydeseq2 not installed. Run: pip install pydeseq2", file=sys.stderr)
        return None

    # Build DataFrames
    count_df = pd.DataFrame(
        count_matrix.T,  # PyDESeq2 wants samples × genes
        index=sample_names,
        columns=gene_names,
    )
    meta_df = pd.DataFrame({"condition": conditions}, index=sample_names)

    # Determine reference level (use the first alphabetically as control)
    unique_conds = sorted(set(conditions))
    if len(unique_conds) < 2:
        return None
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

        stat = DeseqStats(dds, contrast=["condition", unique_conds[1], ref_level], n_cpus=1)
        stat.summary()

        results = stat.results_df
        gene_results = {}
        for gene in results.index:
            row = results.loc[gene]
            gene_results[gene] = {
                "log2fc": float(row["log2FoldChange"]) if not np.isnan(row["log2FoldChange"]) else 0.0,
                "stat": float(row["stat"]) if not np.isnan(row["stat"]) else 0.0,
                "pvalue": float(row["pvalue"]) if not np.isnan(row["pvalue"]) else 1.0,
                "padj": float(row["padj"]) if not np.isnan(row["padj"]) else 1.0,
            }
        return gene_results
    except Exception as e:
        print(f"  WARN: DESeq2 failed: {e}", file=sys.stderr)
        return None


# ============================================================
# GSEA (GSEApy)
# ============================================================

def run_gsea_prerank(gene_results, gene_set_db=GENE_SET_DB):
    """
    Run GSEA prerank on DE results.
    Returns: dict pathway -> NES
    """
    try:
        import gseapy
        import pandas as pd
    except ImportError:
        print("ERROR: gseapy not installed. Run: pip install gseapy", file=sys.stderr)
        return None

    if not gene_results:
        return None

    # Build ranked gene list (use the Wald statistic as ranking metric)
    ranked = {g: r["stat"] for g, r in gene_results.items() if r["stat"] != 0.0}
    if len(ranked) < 50:
        return None

    rnk = pd.Series(ranked).sort_values(ascending=False)

    try:
        import os
        local_gmt = f"/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/{gene_set_db}.gmt"
        actual_gene_set = local_gmt if os.path.exists(local_gmt) else gene_set_db
        pre_res = gseapy.prerank(
            rnk=rnk,
            gene_sets=actual_gene_set,
            outdir=None,         # Don't write files
            min_size=15,
            max_size=500,
            permutation_num=100, # Reduced for speed
            seed=42,
            verbose=False,
            no_plot=True,
        )
        pathway_nes = {}
        for term in pre_res.res2d.index:
            row = pre_res.res2d.loc[term]
            pathway_nes[row["Term"]] = float(row["NES"])
        return pathway_nes
    except Exception as e:
        print(f"  WARN: GSEA failed: {e}", file=sys.stderr)
        return None


# ============================================================
# Concordance Metrics (Section 6)
# ============================================================

def spearman_correlation(vec_a, vec_b):
    """Spearman rank correlation between two aligned vectors."""
    from scipy.stats import spearmanr
    valid = [(a, b) for a, b in zip(vec_a, vec_b) if not (np.isnan(a) or np.isnan(b))]
    if len(valid) < 10:
        return np.nan
    a, b = zip(*valid)
    rho, _ = spearmanr(a, b)
    return rho


def jaccard(set_a, set_b):
    """Jaccard index between two sets."""
    if not set_a and not set_b:
        return 1.0
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)


def compute_concordance(candidate_de, candidate_gsea, ref_gene_z, ref_path_nes,
                        ref_deg_set, ref_path_set, top_path_signs):
    """
    Compute all concordance metrics between candidate result and reference.
    Returns dict of metric -> value.
    """
    metrics = {}

    # (a) Gene-level rank correlation (Section 6.1a)
    common_genes = sorted(set(candidate_de.keys()) & set(ref_gene_z.keys()))
    if common_genes:
        cand_z = [candidate_de[g]["stat"] for g in common_genes]
        ref_z = [ref_gene_z[g] for g in common_genes]
        metrics["rho_gene"] = spearman_correlation(cand_z, ref_z)
    else:
        metrics["rho_gene"] = np.nan

    # (b) DEG overlap (Section 6.1b)
    cand_degs = set()
    for g, r in candidate_de.items():
        if r["padj"] < DEG_FDR and abs(r["log2fc"]) > DEG_LFC:
            cand_degs.add(g)
    # Removed TOP_K_DEG truncation so Jaccard compares true DEG sets
    metrics["jaccard_deg"] = jaccard(ref_deg_set, cand_degs)

    # (c) Pathway-level rank correlation (Section 6.2c)
    if candidate_gsea and ref_path_nes:
        common_paths = sorted(set(candidate_gsea.keys()) & set(ref_path_nes.keys()))
        if common_paths:
            cand_nes = [candidate_gsea[p] for p in common_paths]
            ref_nes_vals = [ref_path_nes[p] for p in common_paths]
            metrics["rho_path"] = spearman_correlation(cand_nes, ref_nes_vals)
        else:
            metrics["rho_path"] = np.nan
    else:
        metrics["rho_path"] = np.nan

    # (d) Significant pathway overlap (Section 6.2d)
    cand_sig_paths = set()
    if candidate_gsea:
        # We don't have FDR from prerank directly in dict, use absolute NES > threshold
        for p, nes in candidate_gsea.items():
            if abs(nes) > 1.0:  # Approximate significance
                cand_sig_paths.add(p)
    metrics["jaccard_path"] = jaccard(ref_path_set, cand_sig_paths)

    # (e) Direction concordance among top pathways (Section 6.2e)
    if candidate_gsea and top_path_signs:
        agree = 0
        total = 0
        for path, ref_sign in top_path_signs.items():
            if path in candidate_gsea:
                total += 1
                cand_sign = np.sign(candidate_gsea[path])
                if cand_sign == ref_sign:
                    agree += 1
        metrics["dir_concordance"] = agree / total if total > 0 else np.nan
    else:
        metrics["dir_concordance"] = np.nan

    return metrics


# ============================================================
# Core: Process One Target SRR
# ============================================================

def process_target_srr(args):
    """
    For one target SRR in one project:
    1. Build leave-one-out references under each stable method
    2. Compute consensus reference
    3. Evaluate each candidate method
    4. Return concordance row
    """
    target_srr, project, project_dir, all_srrs, condition_map = args

    try:
        other_srrs = [s for s in all_srrs if s != target_srr]
        other_conditions = [condition_map[s] for s in other_srrs]
        target_condition = condition_map[target_srr]

        # Check we have enough samples per group after removal
        from collections import Counter
        group_counts = Counter(other_conditions)
        if any(c < MIN_SAMPLES_PER_GROUP for c in group_counts.values()):
            # Fallback: use full project (Section 12)
            other_srrs = all_srrs
            other_conditions = [condition_map[s] for s in other_srrs]

        # Discover gene universe from untrimmed files
        all_genes = None
        for srr in all_srrs:
            fp = os.path.join(project_dir, MODE_PATTERNS["U"].format(srr=srr))
            if os.path.exists(fp):
                c = load_featurecounts(fp)
                if c:
                    if all_genes is None:
                        all_genes = set(c.keys())
                    else:
                        all_genes |= set(c.keys())
        if all_genes is None:
            return None
        all_genes = sorted(all_genes)

        # ---- Step 1: Leave-one-out reference (Section 4.3) ----
        ref_gene_z_by_method = {}   # method -> {gene: z-score}
        ref_path_nes_by_method = {} # method -> {pathway: NES}

        for method in REF_METHODS:
            gene_names, matrix, valid_srrs = build_count_matrix(
                other_srrs, project_dir, method, all_genes
            )
            if matrix is None or len(valid_srrs) < 4:
                continue
            valid_conditions = [condition_map[s] for s in valid_srrs]
            if len(set(valid_conditions)) < 2:
                continue

            de_results = run_deseq2(gene_names, matrix, valid_srrs, valid_conditions)
            if de_results is None:
                continue
            ref_gene_z_by_method[method] = {g: r["stat"] for g, r in de_results.items()}

            gsea_results = run_gsea_prerank(de_results)
            if gsea_results:
                ref_path_nes_by_method[method] = gsea_results

        if not ref_gene_z_by_method:
            return None

        # ---- Step 2: Consensus reference (Section 4.4) ----
        ref_gene_z = {}
        for gene in all_genes:
            vals = [ref_gene_z_by_method[m].get(gene, np.nan)
                    for m in ref_gene_z_by_method]
            vals = [v for v in vals if not np.isnan(v)]
            ref_gene_z[gene] = np.median(vals) if vals else 0.0

        ref_path_nes = {}
        all_paths = set()
        for m in ref_path_nes_by_method.values():
            all_paths |= set(m.keys())
        for path in all_paths:
            vals = [ref_path_nes_by_method[m].get(path, np.nan)
                    for m in ref_path_nes_by_method]
            vals = [v for v in vals if not np.isnan(v)]
            ref_path_nes[path] = np.median(vals) if vals else 0.0

        # Stable DEGs from reference: genes where sign agrees in >=75% of ALL methods
        ref_deg_set = set()
        for gene in all_genes:
            signs = []
            for m in ref_gene_z_by_method:
                z = ref_gene_z_by_method[m].get(gene, 0.0)
                if abs(z) > 1.96:  # Approximately significant
                    signs.append(np.sign(z))
            if signs and signs.count(signs[0]) / len(ref_gene_z_by_method) >= 0.75:
                ref_deg_set.add(gene)

        # Stable pathways from reference
        ref_path_set = set()
        for path in all_paths:
            signs = []
            for m in ref_path_nes_by_method:
                nes = ref_path_nes_by_method[m].get(path, 0.0)
                if abs(nes) > 1.0:
                    signs.append(np.sign(nes))
            if signs and signs.count(signs[0]) / len(ref_path_nes_by_method) >= 0.75:
                ref_path_set.add(path)

        # Top L reference pathways and their signs
        top_paths = sorted(ref_path_nes.items(), key=lambda x: abs(x[1]), reverse=True)[:TOP_L_PATH]
        top_path_signs = {p: np.sign(nes) for p, nes in top_paths if nes != 0}

        # ---- Step 3: Evaluate each candidate method (Section 5) ----
        baseline = BASELINE_METHOD if BASELINE_METHOD != "auto" else choose_project_baseline(project_dir, other_srrs, REF_METHODS)
        results_row = {"SRR_ID": target_srr, "project_id": project, "baseline_method": baseline}

        for method in ALL_METHODS:
            # Build evaluation matrix: other samples at baseline, target at candidate method
            eval_srrs = other_srrs + [target_srr]
            eval_conditions = [condition_map[s] for s in other_srrs] + [target_condition]

            # Load other samples at baseline
            gene_names_b, matrix_b, valid_b = build_count_matrix(
                other_srrs, project_dir, baseline, all_genes
            )
            if matrix_b is None:
                continue

            # Load target sample at candidate method
            target_pattern = MODE_PATTERNS[method].format(srr=target_srr)
            target_fp = os.path.join(project_dir, target_pattern)
            if not os.path.exists(target_fp):
                for m_col in ["rho_gene", "rho_path", "jaccard_deg", "jaccard_path", "dir_concordance"]:
                    results_row[f"{method}_{m_col}"] = "NA"
                continue
            target_counts = load_featurecounts(target_fp)
            if target_counts is None:
                for m_col in ["rho_gene", "rho_path", "jaccard_deg", "jaccard_path", "dir_concordance"]:
                    results_row[f"{method}_{m_col}"] = "NA"
                continue

            # Add target sample to matrix
            target_vec = np.array([target_counts.get(g, 0) for g in all_genes], dtype=np.int64)
            eval_matrix = np.column_stack([matrix_b, target_vec.reshape(-1, 1)])
            valid_eval_srrs = valid_b + [target_srr]
            eval_conds = [condition_map[s] for s in valid_b] + [target_condition]

            if len(set(eval_conds)) < 2:
                for m_col in ["rho_gene", "rho_path", "jaccard_deg", "jaccard_path", "dir_concordance"]:
                    results_row[f"{method}_{m_col}"] = "NA"
                continue

            # Run DE + GSEA
            de_results = run_deseq2(all_genes, eval_matrix, valid_eval_srrs, eval_conds)
            if de_results is None:
                for m_col in ["rho_gene", "rho_path", "jaccard_deg", "jaccard_path", "dir_concordance"]:
                    results_row[f"{method}_{m_col}"] = "NA"
                continue
                
            deseq2_csv = os.path.join(LOSO_RAW_DIR, f"{project}_{target_srr}_{method}_deseq2.csv")
            with open(deseq2_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "log2fc", "padj", "stat"])
                for g, r in de_results.items():
                    writer.writerow([g, r["log2fc"], r["padj"], r["stat"]])

            gsea_results = run_gsea_prerank(de_results)
            if gsea_results is not None:
                gsea_csv = os.path.join(LOSO_RAW_DIR, f"{project}_{target_srr}_{method}_gsea.csv")
                with open(gsea_csv, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["pathway", "nes"])
                    for p, nes in gsea_results.items():
                        writer.writerow([p, nes])

            # Compute concordance
            concordance = compute_concordance(
                de_results, gsea_results, ref_gene_z, ref_path_nes,
                ref_deg_set, ref_path_set, top_path_signs
            )

            for m_col, val in concordance.items():
                if isinstance(val, float) and np.isnan(val):
                    results_row[f"{method}_{m_col}"] = "NA"
                else:
                    results_row[f"{method}_{m_col}"] = f"{val:.6f}"

        return results_row

    except Exception as e:
        print(f"  ERROR processing {target_srr} in {project}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return None


# ============================================================
# Main
# ============================================================

def load_metadata(metadata_dir, project):
    """Load Run,condition CSV for a project."""
    filepath = os.path.join(metadata_dir, f"{project}.csv")
    if not os.path.exists(filepath):
        return None
    condition_map = {}
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            srr = row["Run"].strip()
            cond = row["condition"].strip()
            condition_map[srr] = cond
    return condition_map


def main():
    global GENE_SET_DB, BASELINE_METHOD, REF_METHODS, ALL_METHODS, MIN_SAMPLES_PER_GROUP

    parser = argparse.ArgumentParser(
        description="Leave-one-out DE+GSEA concordance evaluation"
    )
    parser.add_argument("--cores", type=int, default=32, help="Parallel cores")
    parser.add_argument("--project", type=str, default="", help="Process only this project")
    parser.add_argument("--metadata-dir", type=str, default=METADATA_DIR,
                        help="Directory containing Run,condition CSV files")
    parser.add_argument("--gene-set", type=str, default=GENE_SET_DB,
                        help="GSEApy gene set database name")
    parser.add_argument("--baseline-method", type=str, default=BASELINE_METHOD,
                        help="Project baseline method for non-target samples: auto,U,A,P5,P10,P20,P35")
    parser.add_argument("--ref-methods", type=str, default=",".join(REF_METHODS),
                        help="Comma-separated reference methods, e.g. U,A,P5,P10")
    parser.add_argument("--all-methods", type=str, default=",".join(ALL_METHODS),
                        help="Comma-separated candidate methods")
    parser.add_argument("--min-samples-per-group", type=int, default=MIN_SAMPLES_PER_GROUP,
                        help="Minimum samples per group after leave-one-out; else use full-project fallback")
    args = parser.parse_args()

    GENE_SET_DB = args.gene_set
    BASELINE_METHOD = args.baseline_method
    REF_METHODS = [m.strip() for m in args.ref_methods.split(",") if m.strip()]
    ALL_METHODS = [m.strip() for m in args.all_methods.split(",") if m.strip()]
    MIN_SAMPLES_PER_GROUP = args.min_samples_per_group

    os.makedirs(CONCORDANCE_DIR, exist_ok=True)

    # Discover projects with metadata
    project_dirs = sorted(glob.glob(os.path.join(COUNTS_BASE, "*/")))
    projects_to_run = []

    for pdir in project_dirs:
        project = os.path.basename(pdir.rstrip("/"))
        if args.project and project != args.project:
            continue
        condition_map = load_metadata(args.metadata_dir, project)
        if condition_map is None:
            continue
        projects_to_run.append((project, pdir, condition_map))

    if not projects_to_run:
        print("No projects with condition metadata found.")
        print(f"Place CSV files (Run,condition) in: {args.metadata_dir}")
        return

    print(f"Found {len(projects_to_run)} projects with condition metadata")

    for project, pdir, condition_map in projects_to_run:
        print(f"\n{'='*60}")
        print(f"Project: {project} ({len(condition_map)} samples)")
        print(f"{'='*60}")

        # Find SRRs that exist in both metadata and on disk
        available_srrs = []
        for srr in condition_map:
            fp = os.path.join(pdir, MODE_PATTERNS["U"].format(srr=srr))
            if os.path.exists(fp):
                available_srrs.append(srr)

        if len(available_srrs) < 4:
            print(f"  Skipping: only {len(available_srrs)} samples available (need ≥4)")
            continue

        # Check existing output to resume
        out_file = os.path.join(CONCORDANCE_DIR, f"{project}_concordance.tsv")
        completed_srrs = set()
        if os.path.exists(out_file):
            with open(out_file) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    completed_srrs.add(row.get("SRR_ID", ""))
        
        available_srrs = [srr for srr in condition_map.keys() if srr not in completed_srrs]
        total_remaining = len(available_srrs)

        print(f"  Target SRRs to process: {len(available_srrs)}")
        if not available_srrs:
            continue

        all_project_srrs = list(condition_map.keys())
        work_items = [
            (srr, project, pdir, all_project_srrs, condition_map)
            for srr in available_srrs
        ]

        # Process (use fewer cores per SRR to avoid memory contention)
        effective_cores = min(args.cores, len(work_items))
        print(f"  Processing {len(work_items)} target SRRs with {effective_cores} cores...")

        valid_count = 0
        if effective_cores > 1:
            with Pool(processes=effective_cores) as pool:
                for result in pool.imap_unordered(process_target_srr, work_items):
                    if result is not None:
                        file_exists = os.path.exists(out_file)
                        all_cols = list(result.keys())
                        with open(out_file, "a") as f:
                            if not file_exists:
                                f.write("\t".join(all_cols) + "\n")
                            f.write("\t".join(str(result.get(c, "NA")) for c in all_cols) + "\n")
                        valid_count += 1
                        print(f"  -> Saved results for target SRR: {result.get('target_srr')} ({valid_count}/{len(work_items)})", flush=True)
        else:
            for w in work_items:
                result = process_target_srr(w)
                if result is not None:
                    file_exists = os.path.exists(out_file)
                    all_cols = list(result.keys())
                    with open(out_file, "a") as f:
                        if not file_exists:
                            f.write("\t".join(all_cols) + "\n")
                        f.write("\t".join(str(result.get(c, "NA")) for c in all_cols) + "\n")
                    valid_count += 1
                    print(f"  -> Saved results for target SRR: {result.get('target_srr')} ({valid_count}/{len(work_items)})", flush=True)

        print(f"  Completed: {valid_count}/{len(work_items)} written to {out_file}")

    # Write combined concordance file
    all_results = []
    for project, _, _ in projects_to_run:
        proj_file = os.path.join(CONCORDANCE_DIR, f"{project}_concordance.tsv")
        if os.path.exists(proj_file):
            with open(proj_file) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    all_results.append(row)

    if all_results:
        combined_file = os.path.join(CONCORDANCE_DIR, "all_concordance.tsv")
        all_cols = list(all_results[0].keys())
        with open(combined_file, "w") as f:
            f.write("\t".join(all_cols) + "\n")
            for row in all_results:
                f.write("\t".join(str(row.get(c, "NA")) for c in all_cols) + "\n")
        print(f"\nCombined concordance: {combined_file} ({len(all_results)} rows)")


if __name__ == "__main__":
    main()
