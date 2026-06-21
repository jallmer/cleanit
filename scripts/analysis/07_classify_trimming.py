#!/usr/bin/env python3
"""
07_classify_trimming.py — classify trimming outcomes with technical penalties.

Implements Sections 7–10 of the sample-specific trimming methodology:
  1. Load concordance metrics from the DE/GSEA step
  2. Add information-loss penalties from flattened trimming/count/alignment outputs
  3. Choose t* using the hierarchical rule
  4. Classify each method relative to untrimmed
  5. Compute benefit scores with optional net-benefit penalty

Writes:
  /scratch/hpc-prf-omiks/ja/analysis/trimming_classification.tsv
  /scratch/hpc-prf-omiks/ja/analysis/trimming_benefit.tsv
  /scratch/hpc-prf-omiks/ja/analysis/trimming_penalties.tsv
"""

from __future__ import annotations

import csv
import gzip
import os
import sys
from collections import Counter

import numpy as np

OUT_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis"
CONCORDANCE_DIR = os.path.join(OUT_DIR, "concordance")
COUNTS_BASE = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/flattened_counts"
TRIM_FILE = os.path.join(OUT_DIR, "trimmomatic_detail.tsv")
BOWTIE_FILE = os.path.join(OUT_DIR, "bowtie2_alignment_stats.tsv")

ALL_METHODS = ["U", "A", "P5", "P10", "P20", "P35"]
AGGRESSIVENESS_ORDER = {m: i for i, m in enumerate(ALL_METHODS)}
METHOD_TO_MODE = {
    "U": "untrimmed",
    "A": "adapter_only",
    "P5": "P5",
    "P10": "P10",
    "P20": "P20",
    "P35": "P35",
}
MODE_PATTERNS = {
    "untrimmed": "untrmd_{srr}_fC.txt.gz",
    "adapter_only": "{srr}_trimmomatic_adapter_fC.txt.gz",
    "P5": "{srr}_trimmomatic_P5_fC.txt.gz",
    "P10": "{srr}_trimmomatic_P10_fC.txt.gz",
    "P20": "{srr}_trimmomatic_P20_fC.txt.gz",
    "P35": "{srr}_trimmomatic_P35_fC.txt.gz",
}

# Methodology thresholds
DELTA_TOL = 0.01
EPSILON = 0.02

# Technical damage thresholds used for helpful/harmful calls
MAX_HELPFUL_READ_LOSS = 0.15
MAX_HELPFUL_GENE_LOSS = 0.10
MAX_HELPFUL_ASSIGNED_LOSS = 0.05
MATERIAL_PATHWAY_OVERLAP_DROP = 0.05

# Penalty weights for optional B_net
LAMBDA_READ = 0.50
LAMBDA_GENE = 0.30
LAMBDA_ASSIGNED = 0.20


def safe_float(v):
    if v is None or v == "" or str(v).strip() == "NA":
        return np.nan
    try:
        return float(str(v).rstrip("%"))
    except ValueError:
        return np.nan


def load_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_summary_metrics(project, srr, mode):
    base = os.path.join(COUNTS_BASE, project)
    stem = MODE_PATTERNS[mode].format(srr=srr)
    summary_paths = [
        os.path.join(base, stem.replace(".txt.gz", ".txt.summary.gz")),
        os.path.join(base, stem.replace(".txt.gz", ".txt.summary")),
    ]
    assigned = None
    total = None
    for path in summary_paths:
        if not os.path.exists(path):
            continue
        opener = gzip.open if path.endswith(".gz") else open
        assigned = 0
        total = 0
        with opener(path, "rt") as fh:
            next(fh, None)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                val = int(parts[1])
                total += val
                if parts[0] == "Assigned":
                    assigned = val
        break
    assigned_frac = (assigned / total) if assigned is not None and total else np.nan
    return assigned_frac


def count_detected_genes(project, srr, mode):
    base = os.path.join(COUNTS_BASE, project)
    stem = MODE_PATTERNS[mode].format(srr=srr)
    path = os.path.join(base, stem)
    if not os.path.exists(path):
        return np.nan
    opener = gzip.open if path.endswith(".gz") else open
    detected = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 7 and int(parts[6]) > 0:
                detected += 1
    return float(detected)


def build_penalty_lookup():
    trim_rows = load_tsv(TRIM_FILE)
    bowtie_rows = load_tsv(BOWTIE_FILE)

    trim_map = {}
    for row in trim_rows:
        trim_map[(row["project_id"], row["SRR_ID"], row["mode"])] = row

    bowtie_map = {}
    for row in bowtie_rows:
        bowtie_map[(row["project_id"], row["srr_id"], row["mode"])] = row

    srrs = set((r["project_id"], r["srr_id"]) for r in bowtie_rows)
    
    penalty_lookup = {}
    for project, srr in srrs:
        for method in ["U", "A", "P5", "P10", "P20", "P35"]:
            trim_key = None
            if method == "A": trim_key = "adapter_only"
            elif method != "U": trim_key = method
            
            trim_metrics = trim_map.get((project, srr, trim_key), {}) if trim_key else {}
            
            bow_key = None
            if method == "U":
                if (project, srr, "untrimmed") in bowtie_map:
                    bow_key = "untrimmed"
                elif (project, srr, f"untrmd_{srr}") in bowtie_map:
                    bow_key = f"untrmd_{srr}"
            else:
                mapping = {"A": "adapter"}
                sfx = mapping.get(method, method)
                bow_key = f"{srr}_trimmomatic_{sfx}"
                
            bow_metrics = bowtie_map.get((project, srr, bow_key), {})
            
            assigned_frac = load_summary_metrics(project, srr, METHOD_TO_MODE[method])
            detected_genes = count_detected_genes(project, srr, METHOD_TO_MODE[method])
            
            input_reads = safe_float(trim_metrics.get("input_reads"))
            surviving = safe_float(trim_metrics.get("surviving"))
            
            read_loss = np.nan
            if method == "U":
                read_loss = 0.0
            elif not np.isnan(input_reads) and input_reads > 0 and not np.isnan(surviving):
                read_loss = max(0.0, 1.0 - (surviving / input_reads))
                
            penalty_lookup[(project, srr, method)] = {
                "assigned_frac": assigned_frac,
                "detected_genes": detected_genes,
                "overall_alignment_rate": safe_float(bow_metrics.get("overall_alignment_rate")),
                "read_loss": read_loss,
            }

    for project, srr in srrs:
        base = penalty_lookup.get((project, srr, "U"), {})
        base_genes = safe_float(base.get("detected_genes"))
        base_assigned = safe_float(base.get("assigned_frac"))
        for method in ["U", "A", "P5", "P10", "P20", "P35"]:
            entry = penalty_lookup[(project, srr, method)]
            curr_genes = safe_float(entry.get("detected_genes"))
            curr_assigned = safe_float(entry.get("assigned_frac"))
            
            gene_loss = np.nan
            if not np.isnan(base_genes) and base_genes > 0 and not np.isnan(curr_genes):
                gene_loss = max(0.0, 1.0 - (curr_genes / base_genes))
                
            assigned_frac_loss = np.nan
            if not np.isnan(base_assigned) and base_assigned > 0 and not np.isnan(curr_assigned):
                assigned_frac_loss = max(0.0, base_assigned - curr_assigned)
                
            entry["gene_loss"] = gene_loss
            entry["assigned_frac_loss"] = assigned_frac_loss
            
    return penalty_lookup


def classification_label(delta_path, jaccard_path_t, jaccard_path_u, penalties):
    read_loss = penalties.get("read_loss", np.nan)
    gene_loss = penalties.get("gene_loss", np.nan)
    assigned_loss = penalties.get("assigned_frac_loss", np.nan)
    excessive_loss = (
        (not np.isnan(read_loss) and read_loss > MAX_HELPFUL_READ_LOSS)
        or (not np.isnan(gene_loss) and gene_loss > MAX_HELPFUL_GENE_LOSS)
        or (not np.isnan(assigned_loss) and assigned_loss > MAX_HELPFUL_ASSIGNED_LOSS)
    )

    overlap_drop = False
    if not np.isnan(jaccard_path_t) and not np.isnan(jaccard_path_u):
        overlap_drop = jaccard_path_t < (jaccard_path_u - MATERIAL_PATHWAY_OVERLAP_DROP)

    if np.isnan(delta_path):
        return "NA"
    if delta_path > EPSILON and not overlap_drop and not excessive_loss:
        return "helpful"
    if delta_path < -EPSILON or overlap_drop or excessive_loss:
        return "harmful"
    return "neutral"


def penalty_score(penalties):
    total = 0.0
    for value, weight in [
        (penalties.get("read_loss", np.nan), LAMBDA_READ),
        (penalties.get("gene_loss", np.nan), LAMBDA_GENE),
        (penalties.get("assigned_frac_loss", np.nan), LAMBDA_ASSIGNED),
    ]:
        if not np.isnan(value):
            total += weight * value
    return total


def find_optimal_method(method_metrics):
    helpful_methods = []
    
    for method in ALL_METHODS:
        if method == "U":
            continue
        classification = method_metrics.get(method, {}).get("classification", "NA")
        if classification == "helpful":
            helpful_methods.append(method)

    if not helpful_methods:
        # Default to U if no trimming method is strictly helpful
        return "U", method_metrics.get("U", {})

    rho_paths = {}
    for method in helpful_methods:
        rho = method_metrics.get(method, {}).get("rho_path", np.nan)
        if not np.isnan(rho):
            rho_paths[method] = rho

    if not rho_paths:
        return "U", method_metrics.get("U", {})

    max_rho = max(rho_paths.values())
    candidates = {m: rho for m, rho in rho_paths.items() if rho >= max_rho - DELTA_TOL}

    if len(candidates) == 1:
        best = next(iter(candidates))
        return best, method_metrics.get(best, {})

    valid_jaccards = [
        method_metrics.get(m, {}).get("jaccard_path", -1)
        for m in candidates
        if not np.isnan(method_metrics.get(m, {}).get("jaccard_path", np.nan))
    ]
    best_jaccard = max(valid_jaccards) if valid_jaccards else -1
    if best_jaccard > -1:
        candidates = {
            m: rho
            for m, rho in candidates.items()
            if method_metrics.get(m, {}).get("jaccard_path", -1) >= best_jaccard - 0.01
        }
    if len(candidates) == 1:
        best = next(iter(candidates))
        return best, method_metrics.get(best, {})

    best = min(candidates.keys(), key=lambda m: AGGRESSIVENESS_ORDER.get(m, 99))
    return best, method_metrics.get(best, {})


def main():
    combined_file = os.path.join(CONCORDANCE_DIR, "all_concordance.tsv")
    if not os.path.exists(combined_file):
        print(f"ERROR: {combined_file} not found. Run 06_run_de_gsea.py first.")
        sys.exit(1)

    penalty_lookup = build_penalty_lookup()

    rows = load_tsv(combined_file)
    print(f"Loaded {len(rows)} concordance records")

    classification_rows = []
    benefit_rows = []
    penalty_rows = []
    metric_names = ["rho_gene", "rho_path", "jaccard_deg", "jaccard_path", "dir_concordance"]

    for row in rows:
        srr = row["SRR_ID"]
        project = row["project_id"]
        baseline_method = row.get("baseline_method", "NA")

        method_metrics = {}
        for method in ALL_METHODS:
            metrics = {}
            for metric_name in metric_names:
                metrics[metric_name] = safe_float(row.get(f"{method}_{metric_name}", "NA"))
            mode = METHOD_TO_MODE[method]
            penalties = penalty_lookup.get((project, srr, mode), {})
            metrics.update({
                "read_loss": safe_float(penalties.get("read_loss")),
                "gene_loss": safe_float(penalties.get("gene_loss")),
                "assigned_frac_loss": safe_float(penalties.get("assigned_frac_loss")),
                "assigned_frac": safe_float(penalties.get("assigned_frac")),
                "detected_genes": safe_float(penalties.get("detected_genes")),
                "overall_alignment_rate": safe_float(penalties.get("overall_alignment_rate")),
            })
            method_metrics[method] = metrics
            
            penalty_rows.append({
                "project_id": project,
                "SRR_ID": srr,
                "baseline_method": baseline_method,
                "method": method,
                "read_loss": f"{metrics['read_loss']:.6f}" if not np.isnan(metrics["read_loss"]) else "NA",
                "gene_loss": f"{metrics['gene_loss']:.6f}" if not np.isnan(metrics["gene_loss"]) else "NA",
                "assigned_frac_loss": f"{metrics['assigned_frac_loss']:.6f}" if not np.isnan(metrics["assigned_frac_loss"]) else "NA",
                "assigned_frac": f"{metrics['assigned_frac']:.6f}" if not np.isnan(metrics["assigned_frac"]) else "NA",
                "detected_genes": str(int(metrics["detected_genes"])) if not np.isnan(metrics["detected_genes"]) else "NA",
                "overall_alignment_rate": f"{metrics['overall_alignment_rate']:.6f}" if not np.isnan(metrics["overall_alignment_rate"]) else "NA",
            })
            
        # Evaluate classifications
        rho_path_u = method_metrics.get("U", {}).get("rho_path", np.nan)
        jaccard_path_u = method_metrics.get("U", {}).get("jaccard_path", np.nan)
        for method in ALL_METHODS:
            rho_path_t = method_metrics[method].get("rho_path", np.nan)
            jaccard_path_t = method_metrics[method].get("jaccard_path", np.nan)
            delta_path = rho_path_t - rho_path_u if not np.isnan(rho_path_t) and not np.isnan(rho_path_u) else np.nan
            method_metrics[method]["delta_path"] = delta_path
            method_metrics[method]["classification"] = classification_label(delta_path, jaccard_path_t, jaccard_path_u, method_metrics[method])

        t_star, t_star_metrics = find_optimal_method(method_metrics)
        valid_rhos = [method_metrics[m].get("rho_path", np.nan) for m in ALL_METHODS]
        valid_rhos = [v for v in valid_rhos if not np.isnan(v)]
        max_rho = max(valid_rhos) if valid_rhos else np.nan
        benefit = max_rho - rho_path_u if not np.isnan(max_rho) and not np.isnan(rho_path_u) else np.nan
        t_star_penalty = penalty_score(t_star_metrics)
        benefit_net = benefit - t_star_penalty if not np.isnan(benefit) else np.nan

        cls_row = {
            "SRR_ID": srr,
            "project_id": project,
            "baseline_method": baseline_method,
            "t_star": t_star,
        }

        for method in ALL_METHODS:
            cls_row[f"{method}_class"] = method_metrics[method]["classification"]
            delta_path = method_metrics[method]["delta_path"]
            cls_row[f"{method}_delta_path"] = f"{delta_path:.6f}" if not np.isnan(delta_path) else "NA"
            cls_row[f"{method}_read_loss"] = (
                f"{method_metrics[method]['read_loss']:.6f}" if not np.isnan(method_metrics[method]["read_loss"]) else "NA"
            )
            cls_row[f"{method}_gene_loss"] = (
                f"{method_metrics[method]['gene_loss']:.6f}" if not np.isnan(method_metrics[method]["gene_loss"]) else "NA"
            )
            cls_row[f"{method}_assigned_frac_loss"] = (
                f"{method_metrics[method]['assigned_frac_loss']:.6f}"
                if not np.isnan(method_metrics[method]["assigned_frac_loss"])
                else "NA"
            )
        classification_rows.append(cls_row)

        benefit_rows.append({
            "SRR_ID": srr,
            "project_id": project,
            "baseline_method": baseline_method,
            "t_star": t_star,
            "benefit_B": f"{benefit:.6f}" if not np.isnan(benefit) else "NA",
            "benefit_B_net": f"{benefit_net:.6f}" if not np.isnan(benefit_net) else "NA",
            "t_star_penalty": f"{t_star_penalty:.6f}",
            "rho_path_U": f"{rho_path_u:.6f}" if not np.isnan(rho_path_u) else "NA",
            "max_rho_path": f"{max_rho:.6f}" if not np.isnan(max_rho) else "NA",
            "t_star_read_loss": (
                f"{t_star_metrics['read_loss']:.6f}" if not np.isnan(t_star_metrics.get("read_loss", np.nan)) else "NA"
            ),
            "t_star_gene_loss": (
                f"{t_star_metrics['gene_loss']:.6f}" if not np.isnan(t_star_metrics.get("gene_loss", np.nan)) else "NA"
            ),
            "t_star_assigned_frac_loss": (
                f"{t_star_metrics['assigned_frac_loss']:.6f}"
                if not np.isnan(t_star_metrics.get("assigned_frac_loss", np.nan))
                else "NA"
            ),
        })

    outputs = [
        (os.path.join(OUT_DIR, "trimming_classification.tsv"), classification_rows),
        (os.path.join(OUT_DIR, "trimming_benefit.tsv"), benefit_rows),
        (os.path.join(OUT_DIR, "trimming_penalties.tsv"), penalty_rows),
    ]
    for path, records in outputs:
        if not records:
            continue
        columns = list(records[0].keys())
        with open(path, "w") as fh:
            fh.write("\t".join(columns) + "\n")
            for record in records:
                fh.write("\t".join(str(record.get(col, "NA")) for col in columns) + "\n")
        print(f"Written: {path}")

    all_classes = [
        row.get(f"{method}_class", "NA")
        for row in classification_rows
        for method in ALL_METHODS
        if method != "U"
    ]
    print("\nClassification Summary (non-U methods):")
    for cls, count in Counter(all_classes).most_common():
        print(f"  {cls}: {count}")

    t_star_counts = Counter(row["t_star"] for row in classification_rows)
    print("\nOptimal Method Distribution:")
    for method, count in sorted(t_star_counts.items(), key=lambda item: AGGRESSIVENESS_ORDER.get(item[0], 99)):
        print(f"  {method}: {count}")


if __name__ == "__main__":
    main()
