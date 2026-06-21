#!/usr/bin/env python3
"""
04_summarize_results.py — Aggregate all analysis outputs and generate plots.

Merges quality stats + count-matrix evaluation + trimmomatic stats into
comprehensive per-project and global summary tables, and generates plots.

Usage:
    python3 scripts/analysis/04_summarize_results.py

Reads from: /scratch/hpc-prf-omiks/ja/analysis/
Writes to:  /scratch/hpc-prf-omiks/ja/analysis/global_summary.tsv
            /scratch/hpc-prf-omiks/ja/analysis/per_project_summary.tsv
            /scratch/hpc-prf-omiks/ja/analysis/plots/
"""

import os
import sys
import csv
from collections import defaultdict

OUT_DIR = "/scratch/hpc-prf-omiks/ja/analysis"
PLOTS_DIR = os.path.join(OUT_DIR, "plots")

# Try importing matplotlib; fall back gracefully if unavailable
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    print("WARNING: matplotlib not available. Skipping plots.")


def load_tsv(filepath):
    """Load a TSV file, return list of dicts."""
    if not os.path.exists(filepath):
        print(f"  WARN: {filepath} not found, skipping.")
        return []
    rows = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def safe_float(v):
    """Convert string to float, returning None for NA/empty."""
    if v is None or v == "" or v == "NA" or v.strip() == "":
        return None
    try:
        return float(v.rstrip("%"))
    except ValueError:
        return None


def compute_stats(values):
    """Return mean, std, min, max for a list of numeric values (ignoring None)."""
    clean = [v for v in values if v is not None]
    if not clean:
        return {"mean": "NA", "std": "NA", "min": "NA", "max": "NA", "n": 0}
    arr = np.array(clean) if HAS_PLOT else clean
    if HAS_PLOT:
        return {
            "mean": f"{arr.mean():.6f}",
            "std": f"{arr.std():.6f}",
            "min": f"{arr.min():.6f}",
            "max": f"{arr.max():.6f}",
            "n": len(clean),
        }
    else:
        m = sum(clean) / len(clean)
        variance = sum((x - m) ** 2 for x in clean) / len(clean)
        return {
            "mean": f"{m:.6f}",
            "std": f"{variance**0.5:.6f}",
            "min": f"{min(clean):.6f}",
            "max": f"{max(clean):.6f}",
            "n": len(clean),
        }


def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)

    # Load evaluation data
    eval_file = os.path.join(OUT_DIR, "per_srr_eval.tsv")
    eval_data = load_tsv(eval_file)
    print(f"Loaded {len(eval_data)} SRR evaluation records")

    # Load quality data
    quality_file = os.path.join(OUT_DIR, "per_srr_quality.tsv")
    quality_data = load_tsv(quality_file)
    quality_map = {}
    for row in quality_data:
        srr = row.get("SRR_ID", "")
        if srr:
            quality_map[srr] = row
    print(f"Loaded {len(quality_data)} quality records")

    # Load trimmomatic data
    trim_file = os.path.join(OUT_DIR, "trimmomatic_detail.tsv")
    trim_data = load_tsv(trim_file)
    print(f"Loaded {len(trim_data)} trimmomatic records")

    # ===== GLOBAL SUMMARY =====
    global_file = os.path.join(OUT_DIR, "global_summary.tsv")

    # Merge quality + eval
    merged_header = [
        "SRR_ID", "project_id", "sequence_depth", "Q_mean", "Q_median",
        "read_length_mean",
        "untrmd_adptrTrmd_jsd", "untrmd_P5Trmd_jsd", "untrmd_P10Trmd_jsd",
        "untrmd_P20Trmd_jsd", "untrmd_P35Trmd_jsd", "P5Trmd_P35Trmd_jsd",
        "untrmd_adptrTrmd_pear", "untrmd_P5Trmd_pear", "untrmd_P10Trmd_pear",
        "untrmd_P20Trmd_pear", "untrmd_P35Trmd_pear", "P5Trmd_P35Trmd_pear",
    ]

    with open(global_file, "w") as f:
        f.write("\t".join(merged_header) + "\n")
        for row in eval_data:
            srr = row.get("SRR_ID", "")
            q = quality_map.get(srr, {})
            merged = [
                srr,
                row.get("project_id", ""),
                q.get("sequence_depth", "NA"),
                q.get("Q_mean", "NA"),
                q.get("Q_median", "NA"),
                q.get("read_length_mean", "NA"),
            ]
            for col in merged_header[6:]:
                merged.append(row.get(col, "NA"))
            f.write("\t".join(merged) + "\n")
    print(f"Written: {global_file}")

    # ===== PER-PROJECT SUMMARY =====
    pear_cols = [
        "untrmd_adptrTrmd_pear", "untrmd_P5Trmd_pear", "untrmd_P10Trmd_pear",
        "untrmd_P20Trmd_pear", "untrmd_P35Trmd_pear", "P5Trmd_P35Trmd_pear",
    ]
    jsd_cols = [
        "untrmd_adptrTrmd_jsd", "untrmd_P5Trmd_jsd", "untrmd_P10Trmd_jsd",
        "untrmd_P20Trmd_jsd", "untrmd_P35Trmd_jsd", "P5Trmd_P35Trmd_jsd",
    ]

    project_data = defaultdict(list)
    for row in eval_data:
        project_data[row.get("project_id", "unknown")].append(row)

    proj_summary_file = os.path.join(OUT_DIR, "per_project_summary.tsv")
    with open(proj_summary_file, "w") as f:
        header = ["project_id", "n_srrs"]
        for col in pear_cols:
            header.extend([f"{col}_mean", f"{col}_std"])
        for col in jsd_cols:
            header.extend([f"{col}_mean", f"{col}_std"])
        f.write("\t".join(header) + "\n")

        for project in sorted(project_data.keys()):
            rows = project_data[project]
            line = [project, str(len(rows))]
            for col in pear_cols:
                vals = [safe_float(r.get(col)) for r in rows]
                stats = compute_stats(vals)
                line.extend([stats["mean"], stats["std"]])
            for col in jsd_cols:
                vals = [safe_float(r.get(col)) for r in rows]
                stats = compute_stats(vals)
                line.extend([stats["mean"], stats["std"]])
            f.write("\t".join(line) + "\n")
    print(f"Written: {proj_summary_file}")

    # ===== PLOTS =====
    if not HAS_PLOT:
        print("Skipping plots (matplotlib not available)")
        return

    # --- Plot 1: Pearson Correlation Boxplots (like boxplot_pearson_all_SRRs.png) ---
    fig, ax = plt.subplots(figsize=(14, 7))

    labels = ["vs Adapter", "vs P5", "vs P10", "vs P20", "vs P35", "P5 vs P35"]
    data_for_plot = []
    for col in pear_cols:
        vals = [safe_float(r.get(col)) for r in eval_data]
        clean = [v for v in vals if v is not None]
        data_for_plot.append(clean)

    bp = ax.boxplot(data_for_plot, labels=labels, patch_artist=True, notch=True,
                     widths=0.6, showfliers=True, flierprops=dict(markersize=3, alpha=0.4))

    colors = ["#4FC3F7", "#81C784", "#AED581", "#FFD54F", "#FF8A65", "#E57373"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel("Pearson Correlation", fontsize=13)
    ax.set_xlabel("Comparison (Untrimmed vs Trimmed Mode)", fontsize=13)
    ax.set_title("Pearson Correlation of Gene Counts Across Trimming Modes\n(All SRRs)", fontsize=15)
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, "boxplot_pearson_all_SRRs.png")
    fig.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Written: {plot_file}")

    # --- Plot 2: JSD Boxplots ---
    fig, ax = plt.subplots(figsize=(14, 7))
    jsd_for_plot = []
    for col in jsd_cols:
        vals = [safe_float(r.get(col)) for r in eval_data]
        clean = [v for v in vals if v is not None]
        jsd_for_plot.append(clean)

    bp = ax.boxplot(jsd_for_plot, labels=labels, patch_artist=True, notch=True,
                     widths=0.6, showfliers=True, flierprops=dict(markersize=3, alpha=0.4))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel("Jensen-Shannon Divergence", fontsize=13)
    ax.set_xlabel("Comparison (Untrimmed vs Trimmed Mode)", fontsize=13)
    ax.set_title("JSD Between Gene Count Distributions Across Trimming Modes\n(All SRRs)", fontsize=15)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, "boxplot_jsd_all_SRRs.png")
    fig.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Written: {plot_file}")

    # --- Plot 3: Per-project Pearson heatmap ---
    projects = sorted(project_data.keys())
    comparison_labels = ["Adapter", "P5", "P10", "P20", "P35", "P5vP35"]
    matrix = np.full((len(projects), len(pear_cols)), np.nan)

    for i, proj in enumerate(projects):
        for j, col in enumerate(pear_cols):
            vals = [safe_float(r.get(col)) for r in project_data[proj]]
            clean = [v for v in vals if v is not None]
            if clean:
                matrix[i, j] = np.mean(clean)

    fig, ax = plt.subplots(figsize=(10, max(6, len(projects) * 0.4)))
    im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0.9, vmax=1.0)
    ax.set_xticks(range(len(comparison_labels)))
    ax.set_xticklabels(comparison_labels, fontsize=10)
    ax.set_yticks(range(len(projects)))
    ax.set_yticklabels(projects, fontsize=8)
    ax.set_title("Mean Pearson Correlation by Project and Trimming Mode", fontsize=13)
    fig.colorbar(im, ax=ax, label="Pearson r", shrink=0.8)
    fig.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, "heatmap_pearson_by_project.png")
    fig.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Written: {plot_file}")

    print("\nAll summaries and plots complete.")


if __name__ == "__main__":
    main()
