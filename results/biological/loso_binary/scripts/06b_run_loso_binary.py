#!/usr/bin/env python3
"""
06b_run_loso_binary.py — Leave-one-sample-out DE + GSEA concordance evaluation.

Adapted for execution from the HPC working directory.
Reads:  deseq2_sample_annotation.tsv
Counts: flattened_counts/<project>/
Writes: analysis/loso_binary/results/

Usage:
    python3 06b_run_loso_binary.py [--cores 8] [--project PRJNA...]
"""

import os
import sys
import gzip
import csv
import argparse
import fcntl
import signal
import time
import warnings
import traceback
from collections import defaultdict, Counter
from multiprocessing import Pool
import numpy as np

warnings.filterwarnings("ignore")

# === Paths (local) ===
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
def infer_base_dir(script_dir):
    """Return the project scratch dir for both source and deployed script locations."""
    parent = os.path.dirname(script_dir)
    grandparent = os.path.dirname(parent)
    deployed_name = os.path.basename(parent)
    if (os.path.basename(script_dir) == "scripts" and
            (deployed_name == "loso_binary" or deployed_name.startswith("loso_binary_gsea"))):
        return os.path.dirname(grandparent)
    return script_dir


BASE_DIR = infer_base_dir(SCRIPT_DIR)
COUNTS_BASE = os.path.join(BASE_DIR, "flattened_counts")
ANNOTATION_TSV = os.path.join(BASE_DIR, "deseq2_sample_annotation.tsv")
OUT_DIR = os.path.join(BASE_DIR, "analysis", "loso_binary", "results")

# Trimming modes and their filename patterns
MODE_PATTERNS = {
    "U":   "untrmd_{srr}_fC.txt.gz",
    "A":   "{srr}_trimmomatic_adapter_fC.txt.gz",
    "P5":  "{srr}_trimmomatic_P5_fC.txt.gz",
    "P10": "{srr}_trimmomatic_P10_fC.txt.gz",
    "P20": "{srr}_trimmomatic_P20_fC.txt.gz",
    "P35": "{srr}_trimmomatic_P35_fC.txt.gz",
}

REF_METHODS = ["U", "A", "P5", "P10"]
ALL_METHODS = ["U", "A", "P5", "P10", "P20", "P35"]
RESULT_COLUMNS = ["project", "SRR_ID", "method", "rho_gene", "jaccard_deg",
                  "rho_pathway", "jaccard_pathway", "dir_concordance"]
SKIP_COLUMNS = ["timestamp", "project", "SRR_ID", "stage", "method", "reason",
                "detail", "n_samples", "n_genes", "elapsed_seconds"]

GENE_SET_DB = "MSigDB_Hallmark_2020"

DEG_FDR = 0.05
DEG_LFC = 1.0
TOP_L_PATH = 20
MIN_SAMPLES_PER_GROUP = 3
MAX_DESEQ2_SECONDS = int(os.environ.get("LOSO_MAX_DESEQ2_SECONDS", "7200"))


# ============================================================
# Annotation Loader
# ============================================================

def load_annotation(tsv_path, project_filter=None):
    """
    Load the binary annotation table.
    Returns: dict[project] -> list of {SRR_ID, condition, ...}
    Only includes active (non-excluded) SRRs with group > 0.
    """
    import pandas as pd
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)

    # Filter: only active SRRs
    df = df[(df["exclude"].isna()) | (df["exclude"] == "")]
    df = df[df["group"].astype(int) > 0]

    if project_filter:
        df = df[df["BioProject"] == project_filter]

    # Group by project
    projects = {}
    for proj in df["BioProject"].unique():
        p = df[df["BioProject"] == proj]
        samples = []
        for _, row in p.iterrows():
            samples.append({
                "SRR_ID": row["SRR_ID"],
                "condition": row["condition"],  # "control" or "treatment"
                "condition_original": row.get("condition_original", ""),
            })
        projects[proj] = samples

    return projects


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
    Build a count matrix (genes × samples) for a list of SRRs.
    Returns: (gene_names, count_array, valid_srrs)
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

    if all_genes is None:
        all_genes = sorted(set().union(*[set(c.keys()) for c in sample_counts]))

    matrix = np.zeros((len(all_genes), len(valid_srrs)), dtype=np.int64)
    for j, counts in enumerate(sample_counts):
        for i, gene in enumerate(all_genes):
            matrix[i, j] = counts.get(gene, 0)

    return all_genes, matrix, valid_srrs


# ============================================================
# DE Analysis (PyDESeq2)
# ============================================================

class Deseq2Timeout(RuntimeError):
    pass


def timeout_handler(signum, frame):
    raise Deseq2Timeout(f"DESeq2 exceeded {MAX_DESEQ2_SECONDS} seconds")


def append_skip_row(skip_file, row):
    if not skip_file:
        return
    os.makedirs(os.path.dirname(skip_file), exist_ok=True)
    with open(skip_file, "a+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0, os.SEEK_END)
        file_exists = f.tell() > 0
        if not file_exists:
            f.write("\t".join(SKIP_COLUMNS) + "\n")
        f.write("\t".join(str(row.get(c, "")) for c in SKIP_COLUMNS) + "\n")
        f.flush()
        os.fsync(f.fileno())
        fcntl.flock(f, fcntl.LOCK_UN)


def record_deseq2_skip(skip_file, context, reason, detail, n_samples, n_genes, elapsed):
    row = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "project": context.get("project", ""),
        "SRR_ID": context.get("SRR_ID", ""),
        "stage": context.get("stage", ""),
        "method": context.get("method", ""),
        "reason": reason,
        "detail": detail,
        "n_samples": n_samples,
        "n_genes": n_genes,
        "elapsed_seconds": f"{elapsed:.1f}",
    }
    append_skip_row(skip_file, row)


def load_skip_index(skip_file):
    target_skips = set()
    method_skips = defaultdict(set)
    if not skip_file or not os.path.exists(skip_file):
        return target_skips, method_skips
    try:
        with open(skip_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                project = row.get("project", "")
                srr = row.get("SRR_ID", "")
                stage = row.get("stage", "")
                method = row.get("method", "")
                if not project or not srr:
                    continue
                if stage in {"target", "reference"}:
                    target_skips.add((project, srr))
                elif stage == "candidate" and method:
                    method_skips[(project, srr)].add(method)
    except Exception:
        pass
    return target_skips, method_skips


def run_deseq2(gene_names, count_matrix, sample_names, conditions,
               context=None, skip_file=None):
    """
    Run PyDESeq2 on a count matrix with binary conditions (control/treatment).
    Returns: (results, status). status is ok, invalid, failed, or timeout.
    """
    context = context or {}
    start = time.monotonic()
    try:
        import pandas as pd
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("ERROR: pydeseq2 not installed. Run: pip install pydeseq2", file=sys.stderr)
        sys.exit(1)

    count_df = pd.DataFrame(
        count_matrix.T,
        index=sample_names,
        columns=gene_names,
    )
    meta_df = pd.DataFrame({"condition": conditions}, index=sample_names)

    unique_conds = sorted(set(conditions))
    if len(unique_conds) < 2:
        return None, "invalid"

    # Binary: always contrast treatment vs control
    ref_level = "control"
    contrast_level = "treatment"
    if ref_level not in unique_conds or contrast_level not in unique_conds:
        # Fallback: use alphabetical
        ref_level = unique_conds[0]
        contrast_level = unique_conds[1]

    try:
        dds = DeseqDataSet(
            counts=count_df,
            metadata=meta_df,
            design_factors="condition",
            refit_cooks=True,
            n_cpus=1,
        )
        old_handler = None
        if MAX_DESEQ2_SECONDS > 0:
            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(MAX_DESEQ2_SECONDS)
        try:
            dds.deseq2()
            stat = DeseqStats(dds, contrast=["condition", contrast_level, ref_level], n_cpus=1)
            stat.summary()
        finally:
            if MAX_DESEQ2_SECONDS > 0:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)

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
        return gene_results, "ok"
    except Deseq2Timeout as e:
        elapsed = time.monotonic() - start
        detail = f"{e}; limit={MAX_DESEQ2_SECONDS}s"
        print(f"  WARN: pathological DESeq2 case skipped: {detail}", file=sys.stderr)
        record_deseq2_skip(skip_file, context, "deseq2_timeout", detail,
                           len(sample_names), len(gene_names), elapsed)
        return None, "timeout"
    except Exception as e:
        elapsed = time.monotonic() - start
        detail = str(e).replace("\t", " ").replace("\n", " ")
        print(f"  WARN: DESeq2 failed: {detail}", file=sys.stderr)
        record_deseq2_skip(skip_file, context, "deseq2_failed", detail,
                           len(sample_names), len(gene_names), elapsed)
        return None, "failed"


# ============================================================
# GSEA (GSEApy)
# ============================================================

def resolve_gene_set(gene_set_db=GENE_SET_DB):
    """Resolve the local GMT shipped with the analysis before falling back to a GSEApy name."""
    if os.path.exists(gene_set_db):
        return normalize_gmt_for_gseapy(gene_set_db)
    local_gmt = os.path.join(BASE_DIR, "analysis", f"{gene_set_db}.gmt")
    if os.path.exists(local_gmt):
        return normalize_gmt_for_gseapy(local_gmt)
    return gene_set_db


def normalize_gmt_for_gseapy(gmt_path):
    """GSEApy rejects GMT rows with an empty description field; fill it locally."""
    normalized = f"{gmt_path}.gseapy.gmt"
    try:
        if (os.path.exists(normalized) and
                os.path.getmtime(normalized) >= os.path.getmtime(gmt_path)):
            return normalized

        tmp = f"{normalized}.{os.getpid()}.tmp"
        changed = False
        with open(gmt_path) as src, open(tmp, "w") as dst:
            for line in src:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[1] == "":
                    parts[1] = "na"
                    changed = True
                dst.write("\t".join(parts) + "\n")
        if changed:
            os.replace(tmp, normalized)
            return normalized
        os.remove(tmp)
    except Exception as e:
        print(f"  WARN: failed to normalize GMT for GSEApy: {e}", file=sys.stderr)
    return gmt_path


def run_gsea_prerank(gene_results, gene_set_db=GENE_SET_DB):
    """Run GSEA prerank. Returns: dict pathway -> NES"""
    try:
        import gseapy
        import pandas as pd
    except ImportError:
        print("ERROR: gseapy not installed. Run: pip install gseapy", file=sys.stderr)
        sys.exit(1)

    if not gene_results:
        return None

    ranked = {g: r["stat"] for g, r in gene_results.items() if r["stat"] != 0.0}
    if len(ranked) < 50:
        return None

    rnk = pd.Series(ranked).sort_values(ascending=False)

    try:
        pre_res = gseapy.prerank(
            rnk=rnk,
            gene_sets=resolve_gene_set(gene_set_db),
            outdir=None,
            min_size=15,
            max_size=500,
            permutation_num=100,
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
# Concordance Metrics
# ============================================================

def spearman_correlation(vec_a, vec_b):
    from scipy.stats import spearmanr
    valid = [(a, b) for a, b in zip(vec_a, vec_b) if not (np.isnan(a) or np.isnan(b))]
    if len(valid) < 10:
        return np.nan
    a, b = zip(*valid)
    rho, _ = spearmanr(a, b)
    return rho


def jaccard(set_a, set_b):
    if not set_a and not set_b:
        return 1.0
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)


def compute_concordance(candidate_de, candidate_gsea, ref_gene_z, ref_path_nes,
                        ref_deg_set, ref_path_set, top_path_signs):
    metrics = {}

    common_genes = sorted(set(candidate_de.keys()) & set(ref_gene_z.keys()))
    if common_genes:
        cand_z = [candidate_de[g]["stat"] for g in common_genes]
        ref_z = [ref_gene_z[g] for g in common_genes]
        metrics["rho_gene"] = spearman_correlation(cand_z, ref_z)
    else:
        metrics["rho_gene"] = np.nan

    cand_degs = set()
    for g, r in candidate_de.items():
        if r["padj"] < DEG_FDR and abs(r["log2fc"]) > DEG_LFC:
            cand_degs.add(g)
    metrics["jaccard_deg"] = jaccard(ref_deg_set, cand_degs)

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

    cand_sig_paths = set()
    if candidate_gsea:
        for p, nes in candidate_gsea.items():
            if abs(nes) > 1.0:
                cand_sig_paths.add(p)
    metrics["jaccard_path"] = jaccard(ref_path_set, cand_sig_paths)

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

def read_checkpoint_rows(checkpoint_file):
    if not checkpoint_file or not os.path.exists(checkpoint_file):
        return []
    try:
        with open(checkpoint_file) as f:
            return list(csv.DictReader(f, delimiter="\t"))
    except Exception:
        return []


def append_checkpoint_row(checkpoint_file, row):
    if not checkpoint_file:
        return
    os.makedirs(os.path.dirname(checkpoint_file), exist_ok=True)
    with open(checkpoint_file, "a+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0, os.SEEK_END)
        file_exists = f.tell() > 0
        if not file_exists:
            f.write("\t".join(RESULT_COLUMNS) + "\n")
        vals = []
        for c in RESULT_COLUMNS:
            v = row.get(c, "NA")
            if isinstance(v, float) and np.isnan(v):
                vals.append("NA")
            elif isinstance(v, float):
                vals.append(f"{v:.6f}")
            else:
                vals.append(str(v))
        f.write("\t".join(vals) + "\n")
        f.flush()
        os.fsync(f.fileno())
        fcntl.flock(f, fcntl.LOCK_UN)


def process_target_srr(args):
    """
    For one target SRR in one project:
    1. Build leave-one-out references under each stable method
    2. Compute consensus reference
    3. Evaluate each candidate method
    4. Return list of concordance rows (one per method)
    """
    if len(args) == 6:
        target_srr, project, project_dir, all_srrs, condition_map, out_dir = args
    else:
        target_srr, project, project_dir, all_srrs, condition_map = args
        out_dir = None

    checkpoint_file = None
    skip_file = None
    if out_dir:
        checkpoint_file = os.path.join(out_dir, "_target_checkpoints",
                                       project, f"{target_srr}.tsv")
        skip_file = os.path.join(out_dir, "pathological_loso_skips.tsv")

    try:
        checkpoint_rows = read_checkpoint_rows(checkpoint_file)
        target_skips, method_skips = load_skip_index(skip_file)
        skipped_methods = method_skips.get((project, target_srr), set())
        if (project, target_srr) in target_skips:
            return []

        completed_methods = {row.get("method") for row in checkpoint_rows}
        if len((completed_methods | skipped_methods) & set(ALL_METHODS)) >= len(ALL_METHODS):
            return [row for row in checkpoint_rows if row.get("method") in ALL_METHODS]

        other_srrs = [s for s in all_srrs if s != target_srr]
        other_conditions = [condition_map[s] for s in other_srrs]
        target_condition = condition_map[target_srr]

        # Check we have enough samples per group after removal
        group_counts = Counter(other_conditions)
        if any(c < MIN_SAMPLES_PER_GROUP for c in group_counts.values()):
            other_srrs = all_srrs
            other_conditions = [condition_map[s] for s in other_srrs]

        # Discover gene universe
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
            return []
        all_genes = sorted(all_genes)

        # ---- Step 1: Leave-one-out reference ----
        ref_gene_z_by_method = {}
        ref_path_nes_by_method = {}

        for method in REF_METHODS:
            gene_names, matrix, valid_srrs = build_count_matrix(
                other_srrs, project_dir, method, all_genes
            )
            if matrix is None or len(valid_srrs) < 4:
                continue
            valid_conditions = [condition_map[s] for s in valid_srrs]
            if len(set(valid_conditions)) < 2:
                continue

            context = {
                "project": project,
                "SRR_ID": target_srr,
                "stage": "reference",
                "method": method,
            }
            de_results, de_status = run_deseq2(
                gene_names, matrix, valid_srrs, valid_conditions,
                context=context, skip_file=skip_file
            )
            if de_status in {"timeout", "failed"}:
                return []
            if de_results is None:
                continue
            ref_gene_z_by_method[method] = {g: r["stat"] for g, r in de_results.items()}

            gsea_results = run_gsea_prerank(de_results)
            if gsea_results:
                ref_path_nes_by_method[method] = gsea_results

        if not ref_gene_z_by_method:
            return []

        # ---- Step 2: Consensus reference ----
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

        ref_deg_set = set()
        for gene in all_genes:
            signs = []
            for m in ref_gene_z_by_method:
                z = ref_gene_z_by_method[m].get(gene, 0.0)
                if abs(z) > 1.96:
                    signs.append(np.sign(z))
            if signs and signs.count(signs[0]) / len(ref_gene_z_by_method) >= 0.75:
                ref_deg_set.add(gene)

        ref_path_set = set()
        for path in all_paths:
            signs = []
            for m in ref_path_nes_by_method:
                nes = ref_path_nes_by_method[m].get(path, 0.0)
                if abs(nes) > 1.0:
                    signs.append(np.sign(nes))
            if signs and signs.count(signs[0]) / len(ref_path_nes_by_method) >= 0.75:
                ref_path_set.add(path)

        top_paths = sorted(ref_path_nes.items(), key=lambda x: abs(x[1]), reverse=True)[:TOP_L_PATH]
        top_path_signs = {p: np.sign(nes) for p, nes in top_paths if nes != 0}

        # ---- Step 3: Evaluate each candidate method ----
        # Use untrimmed as baseline for non-target samples (simplification: consistent baseline)
        baseline = "U"
        rows = [row for row in checkpoint_rows if row.get("method") in ALL_METHODS]

        for method in ALL_METHODS:
            if method in completed_methods or method in skipped_methods:
                continue

            gene_names_b, matrix_b, valid_b = build_count_matrix(
                other_srrs, project_dir, baseline, all_genes
            )
            if matrix_b is None:
                continue

            target_pattern = MODE_PATTERNS[method].format(srr=target_srr)
            target_fp = os.path.join(project_dir, target_pattern)
            if not os.path.exists(target_fp):
                continue
            target_counts = load_featurecounts(target_fp)
            if target_counts is None:
                continue

            target_vec = np.array([target_counts.get(g, 0) for g in all_genes], dtype=np.int64)
            eval_matrix = np.column_stack([matrix_b, target_vec.reshape(-1, 1)])
            valid_eval_srrs = valid_b + [target_srr]
            eval_conds = [condition_map[s] for s in valid_b] + [target_condition]

            if len(set(eval_conds)) < 2:
                continue

            context = {
                "project": project,
                "SRR_ID": target_srr,
                "stage": "candidate",
                "method": method,
            }
            de_results, de_status = run_deseq2(
                all_genes, eval_matrix, valid_eval_srrs, eval_conds,
                context=context, skip_file=skip_file
            )
            if de_results is None:
                continue

            gsea_results = run_gsea_prerank(de_results)

            concordance = compute_concordance(
                de_results, gsea_results, ref_gene_z, ref_path_nes,
                ref_deg_set, ref_path_set, top_path_signs
            )

            row = {
                "project": project,
                "SRR_ID": target_srr,
                "method": method,
                "rho_gene": concordance.get("rho_gene", np.nan),
                "jaccard_deg": concordance.get("jaccard_deg", np.nan),
                "rho_pathway": concordance.get("rho_path", np.nan),
                "jaccard_pathway": concordance.get("jaccard_path", np.nan),
                "dir_concordance": concordance.get("dir_concordance", np.nan),
            }
            append_checkpoint_row(checkpoint_file, row)
            rows.append(row)

        return rows

    except Exception as e:
        print(f"  ERROR processing {target_srr} in {project}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return []


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Leave-one-out DE+GSEA concordance (binary annotation)"
    )
    parser.add_argument("--cores", type=int, default=8, help="Parallel cores")
    parser.add_argument("--project", type=str, default="", help="Process only this project")
    parser.add_argument("--target-srr", type=str, default="",
                        help="Process only this target SRR within --project")
    parser.add_argument("--annotation", type=str, default=ANNOTATION_TSV,
                        help="Path to the binary annotation TSV")
    parser.add_argument("--counts-dir", type=str, default=COUNTS_BASE,
                        help="Base directory for flattened counts")
    parser.add_argument("--out-dir", type=str, default=OUT_DIR,
                        help="Output directory for concordance results")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Load annotation
    project_filter = args.project if args.project else None
    projects_data = load_annotation(args.annotation, project_filter)

    if not projects_data:
        print("No projects found in annotation table.")
        return

    print(f"Loaded {len(projects_data)} projects from annotation table")
    print(f"Counts directory: {args.counts_dir}")
    print(f"Output directory: {args.out_dir}")

    combined_file = os.path.join(args.out_dir, "bio_concordance_binary.tsv")
    skip_file = os.path.join(args.out_dir, "pathological_loso_skips.tsv")
    header_cols = RESULT_COLUMNS

    total_processed = 0

    for project in sorted(projects_data.keys()):
        samples = projects_data[project]
        condition_map = {s["SRR_ID"]: s["condition"] for s in samples}
        all_srrs = list(condition_map.keys())
        project_dir = os.path.join(args.counts_dir, project)

        if not os.path.isdir(project_dir):
            print(f"  SKIP {project}: counts directory not found at {project_dir}")
            continue

        # Filter to SRRs with actual count files
        available_srrs = [srr for srr in all_srrs
                          if os.path.exists(os.path.join(project_dir,
                                                          MODE_PATTERNS["U"].format(srr=srr)))]

        if len(available_srrs) < 4:
            print(f"  SKIP {project}: only {len(available_srrs)} SRRs with local counts (need ≥4)")
            continue

        # Use only available SRRs
        condition_map_avail = {s: condition_map[s] for s in available_srrs}
        group_counts = Counter(condition_map_avail.values())

        if len(group_counts) < 2:
            print(f"  SKIP {project}: only {len(group_counts)} condition(s) after filtering")
            continue

        # Check existing output to resume
        proj_out = os.path.join(args.out_dir, f"{project}_concordance.tsv")
        completed_srrs = set()
        skipped_targets, skipped_methods = load_skip_index(skip_file)
        if os.path.exists(proj_out):
            import pandas as pd
            try:
                existing = pd.read_csv(proj_out, sep="\t")
                if {"SRR_ID", "method"}.issubset(existing.columns):
                    method_counts = existing.groupby("SRR_ID")["method"].nunique()
                    for srr, n_methods in method_counts.items():
                        n_skipped = len(skipped_methods.get((project, srr), set()))
                        if n_methods + n_skipped >= len(ALL_METHODS):
                            completed_srrs.add(srr)
                else:
                    completed_srrs = set(existing["SRR_ID"].unique())
            except:
                pass
        completed_srrs |= {srr for skipped_project, srr in skipped_targets
                           if skipped_project == project}
        for srr in available_srrs:
            if len(skipped_methods.get((project, srr), set())) >= len(ALL_METHODS):
                completed_srrs.add(srr)

        remaining_srrs = [s for s in available_srrs if s not in completed_srrs]
        if args.target_srr:
            if args.target_srr not in available_srrs:
                print(f"  SKIP {project}: target SRR {args.target_srr} has no local counts")
                continue
            remaining_srrs = [s for s in remaining_srrs if s == args.target_srr]

        print(f"\n{'='*60}")
        print(f"Project: {project} ({len(available_srrs)} SRRs, "
              f"{group_counts.get('control', 0)}c/{group_counts.get('treatment', 0)}t, "
              f"{len(remaining_srrs)} remaining)")
        print(f"{'='*60}")

        if not remaining_srrs:
            print(f"  Already complete, skipping.")
            continue

        work_items = [
            (srr, project, project_dir, available_srrs, condition_map_avail, args.out_dir)
            for srr in remaining_srrs
        ]

        effective_cores = min(args.cores, len(work_items))
        print(f"  Processing {len(work_items)} target SRRs with {effective_cores} cores...")

        valid_count = 0

        def write_rows(rows, proj_out_path):
            nonlocal valid_count
            if not rows:
                return
            os.makedirs(os.path.dirname(proj_out_path), exist_ok=True)
            with open(proj_out_path, "a+") as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                f.seek(0)
                existing_keys = set()
                try:
                    reader = csv.DictReader(f, delimiter="\t")
                    existing_keys = {(r.get("SRR_ID"), r.get("method")) for r in reader}
                except Exception:
                    existing_keys = set()
                rows_to_write = [
                    row for row in rows
                    if (str(row.get("SRR_ID")), str(row.get("method"))) not in existing_keys
                ]
                if not rows_to_write:
                    fcntl.flock(f, fcntl.LOCK_UN)
                    return
                f.seek(0, os.SEEK_END)
                file_exists = f.tell() > 0
                if not file_exists:
                    f.write("\t".join(header_cols) + "\n")
                for row in rows_to_write:
                    vals = []
                    for c in header_cols:
                        v = row.get(c, "NA")
                        if isinstance(v, float) and np.isnan(v):
                            vals.append("NA")
                        elif isinstance(v, float):
                            vals.append(f"{v:.6f}")
                        else:
                            vals.append(str(v))
                    f.write("\t".join(vals) + "\n")
                f.flush()
                os.fsync(f.fileno())
                fcntl.flock(f, fcntl.LOCK_UN)
            valid_count += len(rows_to_write)

        if effective_cores > 1:
            with Pool(processes=effective_cores) as pool:
                for result_rows in pool.imap_unordered(process_target_srr, work_items):
                    if result_rows:
                        write_rows(result_rows, proj_out)
                        srr_id = result_rows[0]["SRR_ID"] if result_rows else "?"
                        n_done = valid_count // len(ALL_METHODS)
                        print(f"  -> {srr_id} done ({n_done}/{len(work_items)})", flush=True)
        else:
            for w in work_items:
                result_rows = process_target_srr(w)
                if result_rows:
                    write_rows(result_rows, proj_out)
                    srr_id = result_rows[0]["SRR_ID"] if result_rows else "?"
                    n_done = valid_count // len(ALL_METHODS)
                    print(f"  -> {srr_id} done ({n_done}/{len(work_items)})", flush=True)

        total_processed += valid_count
        print(f"  Completed: {valid_count} rows written to {proj_out}")

    if args.target_srr:
        print("\nSingle-target run complete; skipping combined file rebuild.")
        print(f"\nTotal rows written: {total_processed}")
        return

    # Combine all per-project files
    print(f"\n{'='*60}")
    print(f"Combining results...")
    all_rows = []
    for project in sorted(projects_data.keys()):
        proj_file = os.path.join(args.out_dir, f"{project}_concordance.tsv")
        if os.path.exists(proj_file):
            with open(proj_file) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    all_rows.append(row)

    if all_rows:
        with open(combined_file, "w") as f:
            f.write("\t".join(header_cols) + "\n")
            for row in all_rows:
                f.write("\t".join(str(row.get(c, "NA")) for c in header_cols) + "\n")
        print(f"Combined: {combined_file} ({len(all_rows)} rows)")
    else:
        print("No results to combine.")

    print(f"\nTotal rows written: {total_processed}")


if __name__ == "__main__":
    main()
