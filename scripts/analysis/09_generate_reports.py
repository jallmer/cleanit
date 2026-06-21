#!/usr/bin/env python3
"""
09_generate_reports.py — Publication-quality figures and summary tables.

Aggregates all analysis outputs into final tables and plots.

Usage:
    python3 scripts/analysis/09_generate_reports.py

Reads from:  /scratch/hpc-prf-omiks/ja/analysis/
Writes to:   /scratch/hpc-prf-omiks/ja/analysis/plots/
             /scratch/hpc-prf-omiks/ja/analysis/final_summary.tsv
"""

import os
import sys
import csv
from collections import defaultdict, Counter
import numpy as np

OUT_DIR = "/scratch/hpc-prf-omiks/ja/analysis"
PLOTS_DIR = os.path.join(OUT_DIR, "plots")
CONCORDANCE_DIR = os.path.join(OUT_DIR, "concordance")

ALL_METHODS = ["U", "A", "P5", "P10", "P20", "P35"]
METHOD_LABELS = {
    "U": "Untrimmed", "A": "Adapter", "P5": "Phred 5",
    "P10": "Phred 10", "P20": "Phred 20", "P35": "Phred 35",
}

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    print("WARNING: matplotlib not available. Skipping plots.")


def safe_float(v):
    if v is None or v == "" or v.strip() == "NA":
        return np.nan
    try:
        return float(v.rstrip("%"))
    except ValueError:
        return np.nan


def load_tsv(filepath):
    if not os.path.exists(filepath):
        return []
    with open(filepath) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)

    # Load all available data
    eval_data = load_tsv(os.path.join(OUT_DIR, "per_srr_eval.tsv"))
    quality_data = load_tsv(os.path.join(OUT_DIR, "per_srr_quality.tsv"))
    classification_data = load_tsv(os.path.join(OUT_DIR, "trimming_classification.tsv"))
    benefit_data = load_tsv(os.path.join(OUT_DIR, "trimming_benefit.tsv"))
    feature_data = load_tsv(os.path.join(OUT_DIR, "sample_feature_table.tsv"))
    concordance_data = load_tsv(os.path.join(CONCORDANCE_DIR, "all_concordance.tsv"))

    print(f"Data loaded:")
    print(f"  Eval (JSD/Pearson):    {len(eval_data)} SRRs")
    print(f"  Quality (FastQC):      {len(quality_data)} SRRs")
    print(f"  Classification:        {len(classification_data)} SRRs")
    print(f"  Benefit scores:        {len(benefit_data)} SRRs")
    print(f"  Concordance (DE/GSEA): {len(concordance_data)} SRRs")
    print(f"  Feature table:         {len(feature_data)} SRRs")

    # ============================================================
    # Report 1: Classification Summary Table
    # ============================================================
    if classification_data:
        print("\n" + "=" * 60)
        print("Classification Summary")
        print("=" * 60)

        # Per-project breakdown
        project_classifications = defaultdict(lambda: defaultdict(lambda: Counter()))
        for row in classification_data:
            proj = row["project_id"]
            for method in ALL_METHODS:
                cls = row.get(f"{method}_class", "NA")
                project_classifications[proj][method][cls] += 1

        summary_file = os.path.join(OUT_DIR, "classification_summary.tsv")
        with open(summary_file, "w") as f:
            header = ["project_id", "n_srrs"] + [
                f"{m}_{c}" for m in ALL_METHODS[1:]  # Skip U
                for c in ["helpful", "neutral", "harmful"]
            ]
            f.write("\t".join(header) + "\n")
            for proj in sorted(project_classifications.keys()):
                n_srrs = sum(project_classifications[proj]["A"].values())
                line = [proj, str(n_srrs)]
                for method in ALL_METHODS[1:]:
                    counts = project_classifications[proj][method]
                    for cls in ["helpful", "neutral", "harmful"]:
                        line.append(str(counts.get(cls, 0)))
                f.write("\t".join(line) + "\n")
        print(f"  Written: {summary_file}")

        # Global counts
        global_cls = Counter()
        for row in classification_data:
            for method in ALL_METHODS[1:]:
                cls = row.get(f"{method}_class", "NA")
                global_cls[cls] += 1
        print(f"  Global: {dict(global_cls)}")

    # ============================================================
    # Report 1b: Benefit summary
    # ============================================================
    if benefit_data:
        ben_file = os.path.join(OUT_DIR, "benefit_summary.tsv")
        per_project = defaultdict(list)
        with open(ben_file, "w") as f:
            header = [
                "project_id", "n_srrs",
                "benefit_B_mean", "benefit_B_median",
                "benefit_B_net_mean", "benefit_B_net_median",
            ]
            f.write("\t".join(header) + "\n")
            for row in benefit_data:
                per_project[row["project_id"]].append(row)
            for proj in sorted(per_project):
                vals = [safe_float(r.get("benefit_B", "NA")) for r in per_project[proj]]
                vals = [v for v in vals if not np.isnan(v)]
                vals_net = [safe_float(r.get("benefit_B_net", "NA")) for r in per_project[proj]]
                vals_net = [v for v in vals_net if not np.isnan(v)]
                line = [
                    proj,
                    str(len(per_project[proj])),
                    f"{np.mean(vals):.4f}" if vals else "NA",
                    f"{np.median(vals):.4f}" if vals else "NA",
                    f"{np.mean(vals_net):.4f}" if vals_net else "NA",
                    f"{np.median(vals_net):.4f}" if vals_net else "NA",
                ]
                f.write("\t".join(line) + "\n")
        print(f"  Written: {ben_file}")

    # ============================================================
    # Report 2: Concordance Summary by Method
    # ============================================================
    if concordance_data:
        metrics = ["rho_gene", "rho_path", "jaccard_deg", "jaccard_path", "dir_concordance"]
        conc_summary_file = os.path.join(OUT_DIR, "concordance_summary_by_method.tsv")
        with open(conc_summary_file, "w") as f:
            header = ["method"] + [f"{m}_{s}" for m in metrics for s in ["mean", "std", "min", "max"]]
            f.write("\t".join(header) + "\n")
            for method in ALL_METHODS:
                line = [method]
                for metric in metrics:
                    key = f"{method}_{metric}"
                    vals = [safe_float(r.get(key, "NA")) for r in concordance_data]
                    clean = [v for v in vals if not np.isnan(v)]
                    if clean:
                        arr = np.array(clean)
                        line.extend([f"{arr.mean():.4f}", f"{arr.std():.4f}",
                                     f"{arr.min():.4f}", f"{arr.max():.4f}"])
                    else:
                        line.extend(["NA", "NA", "NA", "NA"])
                f.write("\t".join(line) + "\n")
        print(f"  Written: {conc_summary_file}")

    # ============================================================
    # Plots
    # ============================================================
    if not HAS_PLOT:
        print("\nSkipping plots (matplotlib not available).")
        return

    # --- Plot 1: Concordance boxplots by method ---
    if concordance_data:
        metrics_to_plot = ["rho_path", "rho_gene", "jaccard_deg", "dir_concordance"]
        metric_labels = {
            "rho_path": "Pathway\nSpearman ρ",
            "rho_gene": "Gene-Level\nSpearman ρ",
            "jaccard_deg": "DEG\nJaccard",
            "dir_concordance": "Direction\nConcordance",
        }

        fig, axes = plt.subplots(1, len(metrics_to_plot), figsize=(20, 6), sharey=False)
        colors = ["#4FC3F7", "#81C784", "#AED581", "#FFD54F", "#FF8A65", "#E57373"]

        for ax_idx, metric in enumerate(metrics_to_plot):
            ax = axes[ax_idx]
            data_for_plot = []
            for method in ALL_METHODS:
                key = f"{method}_{metric}"
                vals = [safe_float(r.get(key, "NA")) for r in concordance_data]
                clean = [v for v in vals if not np.isnan(v)]
                data_for_plot.append(clean)

            if any(len(d) > 0 for d in data_for_plot):
                bp = ax.boxplot(
                    data_for_plot,
                    labels=[METHOD_LABELS[m] for m in ALL_METHODS],
                    patch_artist=True, notch=True, widths=0.6,
                    showfliers=True,
                    flierprops=dict(markersize=3, alpha=0.4),
                )
                for patch, color in zip(bp["boxes"], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

            ax.set_title(metric_labels[metric], fontsize=11)
            ax.tick_params(axis="x", rotation=45)
            ax.grid(axis="y", alpha=0.3)

        fig.suptitle("Concordance Metrics Across Trimming Methods", fontsize=14, y=1.02)
        fig.tight_layout()
        plot_file = os.path.join(PLOTS_DIR, "concordance_boxplots.png")
        fig.savefig(plot_file, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Written: {plot_file}")

    # --- Plot 2: Classification stacked bar ---
    if classification_data:
        projects = sorted(set(r["project_id"] for r in classification_data))
        cls_colors = {"helpful": "#66BB6A", "neutral": "#FFA726", "harmful": "#EF5350", "NA": "#BDBDBD"}
        fig, axes = plt.subplots(1, len(ALL_METHODS) - 1, figsize=(22, 6), sharey=True)

        for ax_idx, method in enumerate(ALL_METHODS[1:]):
            ax = axes[ax_idx]
            bottoms = np.zeros(len(projects))
            for cls in ["helpful", "neutral", "harmful"]:
                counts = []
                for proj in projects:
                    cnt = sum(1 for r in classification_data
                              if r["project_id"] == proj and r.get(f"{method}_class") == cls)
                    counts.append(cnt)
                ax.bar(range(len(projects)), counts, bottom=bottoms,
                       color=cls_colors[cls], label=cls if ax_idx == 0 else "",
                       edgecolor="white", linewidth=0.5)
                bottoms += np.array(counts)

            ax.set_title(METHOD_LABELS[method], fontsize=11)
            ax.set_xticks(range(len(projects)))
            ax.set_xticklabels(projects, rotation=90, fontsize=7)
            if ax_idx == 0:
                ax.set_ylabel("Number of SRRs", fontsize=11)

        handles = [Patch(facecolor=cls_colors[c], label=c.capitalize()) for c in ["helpful", "neutral", "harmful"]]
        fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=10, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle("Trimming Classification by Project and Method", fontsize=14, y=1.06)
        fig.tight_layout()
        plot_file = os.path.join(PLOTS_DIR, "classification_stacked_bar.png")
        fig.savefig(plot_file, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Written: {plot_file}")

    # --- Plot 2b: Benefit distributions ---
    if benefit_data:
        methods = ["U", "A", "P5", "P10", "P20", "P35"]
        method_colors = {
            "U": "#4FC3F7", "A": "#81C784", "P5": "#AED581",
            "P10": "#FFD54F", "P20": "#FF8A65", "P35": "#E57373"
        }
        t_star_counts = Counter(r.get("t_star", "NA") for r in benefit_data)
        raw_vals = [safe_float(r.get("benefit_B", "NA")) for r in benefit_data]
        raw_vals = [v for v in raw_vals if not np.isnan(v)]
        net_vals = [safe_float(r.get("benefit_B_net", "NA")) for r in benefit_data]
        net_vals = [v for v in net_vals if not np.isnan(v)]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        ax1.hist(raw_vals, bins=25, color="#5C8D89", alpha=0.85, label="benefit_B")
        if net_vals:
            ax1.hist(net_vals, bins=25, color="#C76D3A", alpha=0.55, label="benefit_B_net")
        ax1.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
        ax1.set_title("Benefit Score Distribution")
        ax1.set_xlabel("Benefit")
        ax1.set_ylabel("Count")
        ax1.legend()
        ax1.grid(alpha=0.3)

        ax2.bar(
            methods,
            [t_star_counts.get(m, 0) for m in methods],
            color=[method_colors[m] for m in methods],
            edgecolor="white",
            linewidth=0.6,
        )
        ax2.set_title("Optimal Method Counts")
        ax2.set_xlabel("t*")
        ax2.set_ylabel("SRRs")
        ax2.grid(axis="y", alpha=0.3)

        fig.tight_layout()
        plot_file = os.path.join(PLOTS_DIR, "benefit_summary.png")
        fig.savefig(plot_file, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Written: {plot_file}")

    # --- Plot 3: Pearson + JSD combined from Tier 1 (if available) ---
    if eval_data:
        pear_cols = [
            "untrmd_adptrTrmd_pear", "untrmd_P5Trmd_pear", "untrmd_P10Trmd_pear",
            "untrmd_P20Trmd_pear", "untrmd_P35Trmd_pear", "P5Trmd_P35Trmd_pear",
        ]
        jsd_cols = [
            "untrmd_adptrTrmd_jsd", "untrmd_P5Trmd_jsd", "untrmd_P10Trmd_jsd",
            "untrmd_P20Trmd_jsd", "untrmd_P35Trmd_jsd", "P5Trmd_P35Trmd_jsd",
        ]
        col_labels = ["vs Adapter", "vs P5", "vs P10", "vs P20", "vs P35", "P5 vs P35"]
        colors = ["#4FC3F7", "#81C784", "#AED581", "#FFD54F", "#FF8A65", "#E57373"]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 7))

        # Pearson
        pear_data = []
        for col in pear_cols:
            vals = [safe_float(r.get(col)) for r in eval_data]
            pear_data.append([v for v in vals if not np.isnan(v)])
        bp1 = ax1.boxplot(pear_data, labels=col_labels, patch_artist=True, notch=True, widths=0.6,
                          flierprops=dict(markersize=3, alpha=0.4))
        for patch, c in zip(bp1["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.7)
        ax1.set_ylabel("Pearson Correlation", fontsize=12)
        ax1.set_title("Pearson Correlation (All SRRs)", fontsize=13)
        ax1.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
        ax1.grid(axis="y", alpha=0.3)

        # JSD
        jsd_data = []
        for col in jsd_cols:
            vals = [safe_float(r.get(col)) for r in eval_data]
            jsd_data.append([v for v in vals if not np.isnan(v)])
        bp2 = ax2.boxplot(jsd_data, labels=col_labels, patch_artist=True, notch=True, widths=0.6,
                          flierprops=dict(markersize=3, alpha=0.4))
        for patch, c in zip(bp2["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.7)
        ax2.set_ylabel("Jensen-Shannon Divergence", fontsize=12)
        ax2.set_title("JSD (All SRRs)", fontsize=13)
        ax2.grid(axis="y", alpha=0.3)

        fig.suptitle("Count-Level Trimming Impact Across All SRRs", fontsize=15, y=1.02)
        fig.tight_layout()
        plot_file = os.path.join(PLOTS_DIR, "combined_pearson_jsd_boxplots.png")
        fig.savefig(plot_file, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Written: {plot_file}")

    # --- Plot 4: Per-project heatmap: Δ_path by method ---
    if classification_data:
        projects = sorted(set(r["project_id"] for r in classification_data))
        delta_methods = ALL_METHODS[1:]  # Skip U (delta is 0 by definition)
        matrix = np.full((len(projects), len(delta_methods)), np.nan)

        project_data = defaultdict(list)
        for row in classification_data:
            project_data[row["project_id"]].append(row)

        for i, proj in enumerate(projects):
            for j, method in enumerate(delta_methods):
                key = f"{method}_delta_path"
                vals = [safe_float(r.get(key, "NA")) for r in project_data[proj]]
                clean = [v for v in vals if not np.isnan(v)]
                if clean:
                    matrix[i, j] = np.mean(clean)

        fig, ax = plt.subplots(figsize=(10, max(6, len(projects) * 0.5)))
        vmax = max(0.05, np.nanmax(np.abs(matrix)))
        im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(delta_methods)))
        ax.set_xticklabels([METHOD_LABELS[m] for m in delta_methods], fontsize=10)
        ax.set_yticks(range(len(projects)))
        ax.set_yticklabels(projects, fontsize=9)
        ax.set_title("Mean Δρ_path (vs Untrimmed) by Project and Method", fontsize=13)
        fig.colorbar(im, ax=ax, label="Δρ_path", shrink=0.8)

        # Annotate cells
        for i in range(len(projects)):
            for j in range(len(delta_methods)):
                val = matrix[i, j]
                if not np.isnan(val):
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center", fontsize=7,
                            color="white" if abs(val) > vmax * 0.6 else "black")

        fig.tight_layout()
        plot_file = os.path.join(PLOTS_DIR, "delta_path_heatmap.png")
        fig.savefig(plot_file, dpi=150)
        plt.close()
        print(f"  Written: {plot_file}")

    print("\nAll reports generated.")


if __name__ == "__main__":
    main()
