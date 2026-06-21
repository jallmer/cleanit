#!/usr/bin/env python3
"""
11_compute_throughput.py — Calculate MB/s throughput for each pipeline stage.

Combines timing data (from flattened_timings/) with FASTQ file sizes to compute
throughput in MB/s for: FastQC, trimming (per PHRED level), alignment, counting.

Usage:
    python3 scripts/analysis/11_compute_throughput.py [--cores 32] [--project PRJNA...]

Reads from:
    /scratch/hpc-prf-omiks/ja/flattened_timings/<PROJECT>/<SRR>_timings.tsv
    /scratch/hpc-prf-omiks/ja/flattened_fastqc_raw/<PROJECT>/  (for FASTQ sizes via FastQC)
    /scratch/hpc-prf-omiks/ja/analysis/per_srr_quality.tsv     (for read count + length)
Writes to:
    /scratch/hpc-prf-omiks/ja/analysis/throughput_detail.tsv
    /scratch/hpc-prf-omiks/ja/analysis/throughput_summary.tsv
    /scratch/hpc-prf-omiks/ja/analysis/plots/throughput_boxplots.png
"""

import os
import sys
import csv
import glob
import argparse
from collections import defaultdict
from multiprocessing import Pool
import gzip as gzip_mod
import sqlite3
import numpy as np

TIMINGS_BASE = "/scratch/hpc-prf-omiks/ja/flattened_timings"
COUNTS_BASE = "/scratch/hpc-prf-omiks/ja/flattened_counts"
FLATTENED_BOWTIE_BASE = "/scratch/hpc-prf-omiks/ja/flattened_bowtie2_stats"
DB_FILE = os.path.expanduser("~/srr_queue.db")
QUALITY_FILE = "/scratch/hpc-prf-omiks/ja/analysis/per_srr_quality.tsv"
OUT_DIR = "/scratch/hpc-prf-omiks/ja/analysis"
PLOTS_DIR = os.path.join(OUT_DIR, "plots")

# Stages we want to compute throughput for
# Timing file stages look like: fastqc, trim_P5, align_P5, count_P5, etc.
STAGE_CATEGORIES = {
    "fastqc":  ["fastqc"],
    "trim":    ["trim_adapter_only", "trim_P5", "trim_P10", "trim_P20", "trim_P35"],
    "align":   ["align_untrimmed", "align_adapter_only", "align_P5", "align_P10", "align_P20", "align_P35"],
    "count":   ["count_untrimmed", "count_adapter_only", "count_P5", "count_P10", "count_P20", "count_P35"],
}

# For per-mode breakdown
TRIM_MODES = ["adapter_only", "P5", "P10", "P20", "P35"]
ALIGN_MODES = ["untrimmed", "adapter_only", "P5", "P10", "P20", "P35"]

# Cores used per stage (from v2_master_node.sh: ACTIVE_THREADS=36, FastQC -t 8)
CORES_PER_STAGE = {
    "fastqc": 8,
    "sra_conversion": 36,
    "trim": 36,
    "align": 36,
    "count": 36,
    "mode_total": 36,
    "pipeline_total": 48,  # Full SLURM allocation
    "other": 1,
}


def safe_float(v):
    if v is None or v == "":
        return np.nan
    try:
        return float(v)
    except ValueError:
        return np.nan


def load_timing_file(filepath):
    """Load a timing TSV for one SRR. Returns dict stage -> duration_sec."""
    stages = {}
    try:
        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                stage = row.get("stage", "")
                dur = safe_float(row.get("duration_sec", ""))
                if stage and not np.isnan(dur):
                    stages[stage] = dur
    except Exception as e:
        return {}
    return stages


def estimate_fastq_mb(total_reads, mean_read_length):
    """
    Estimate uncompressed FASTQ file size in MB from read count and length.
    Each read in FASTQ has 4 lines: @header, sequence, +, quality.
    Approximate bytes per read ≈ (header~40) + read_length + 2 + read_length + newlines ≈ 2*read_length + 50
    Compressed (gzip) is typically ~25-30% of uncompressed for FASTQ.
    We report UNCOMPRESSED throughput to make numbers meaningful.
    """
    if np.isnan(total_reads) or np.isnan(mean_read_length):
        return np.nan
    bytes_per_read = 2 * mean_read_length + 50  # sequence + quality + headers + newlines
    total_bytes = total_reads * bytes_per_read
    return total_bytes / (1024 * 1024)  # MB


def get_total_reads_from_summary(project_dir, srr):
    """
    Extract total reads from the untrimmed featureCounts summary file.
    The summary lists: Status<tab>count for categories like
    Assigned, Unassigned_Unmapped, etc. Total = sum of all counts.
    """
    patterns = [
        os.path.join(project_dir, f"untrmd_{srr}_fC.txt.summary.gz"),
        os.path.join(project_dir, f"untrmd_{srr}_fC.txt.summary"),
    ]
    for summary_path in patterns:
        if not os.path.exists(summary_path):
            continue
        try:
            opener = gzip_mod.open if summary_path.endswith(".gz") else open
            total = 0
            with opener(summary_path, "rt") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2 and parts[0] != "Status":
                        try:
                            total += int(parts[1])
                        except ValueError:
                            pass
            if total > 0:
                return total
        except Exception:
            pass
    return None


def load_db_sizes():
    """
    Load size_gb from the SQLite DB as a rough proxy for data volume.
    Returns dict srr_id -> size_gb.
    """
    sizes = {}
    if not os.path.exists(DB_FILE):
        return sizes
    try:
        conn = sqlite3.connect(DB_FILE)
        cursor = conn.execute("SELECT srr_id, size_gb FROM srr_queue WHERE size_gb IS NOT NULL")
        for row in cursor:
            sizes[row[0]] = row[1]
        conn.close()
    except Exception:
        pass
    return sizes


def build_quality_map_from_counts():
    """
    Build a quality map by scanning featureCounts summaries for total reads.
    Assumes a default read length of 150 bp when unknown.
    """
    DEFAULT_READ_LENGTH = 150
    quality_map = {}
    if not os.path.isdir(COUNTS_BASE):
        return quality_map

    for project_dir in sorted(glob.glob(os.path.join(COUNTS_BASE, "*/"))):
        project = os.path.basename(project_dir.rstrip("/"))
        # Find untrimmed summary files to get total reads
        for summary_file in glob.glob(os.path.join(project_dir, "untrmd_*_fC.txt.summary*")):
            basename = os.path.basename(summary_file)
            srr = basename.replace("untrmd_", "").replace("_fC.txt.summary.gz", "").replace("_fC.txt.summary", "")
            if srr in quality_map:
                continue
            total_reads = get_total_reads_from_summary(project_dir, srr)
            if total_reads and total_reads > 0:
                quality_map[srr] = {
                    "SRR_ID": srr,
                    "project_id": project,
                    "sequence_depth": str(total_reads),
                    "read_length_mean": str(DEFAULT_READ_LENGTH),  # Assumed default
                }

    return quality_map


def supplement_quality_from_flattened_bowtie(quality_map):
    added = 0
    if not os.path.isdir(FLATTENED_BOWTIE_BASE):
        return added
    for bowtie_file in glob.glob(os.path.join(FLATTENED_BOWTIE_BASE, "*", "*_bowtie2_stats.tsv")):
        try:
            with open(bowtie_file) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    srr = row.get("srr_id", "")
                    mode = row.get("mode", "")
                    total = safe_float(row.get("total_reads", ""))
                    if mode != "untrimmed":
                        continue
                    if srr and not np.isnan(total) and srr not in quality_map:
                        quality_map[srr] = {
                            "SRR_ID": srr,
                            "project_id": row.get("project_id", ""),
                            "sequence_depth": str(int(total)),
                            "read_length_mean": "150",
                        }
                        added += 1
                        break
        except Exception:
            pass
    return added


def process_project(args):
    """Process all SRRs in a project directory."""
    project, timing_dir, quality_map = args
    results = []

    timing_files = glob.glob(os.path.join(timing_dir, "*_timings.tsv"))

    for tf in timing_files:
        basename = os.path.basename(tf)
        srr = basename.replace("_timings.tsv", "")

        stages = load_timing_file(tf)
        if not stages:
            continue

        # Get quality info for size estimation
        q = quality_map.get(srr, {})
        total_reads = safe_float(q.get("sequence_depth", ""))
        mean_rl = safe_float(q.get("read_length_mean", ""))
        fastq_mb = estimate_fastq_mb(total_reads, mean_rl)

        # Compute throughput for each stage
        for stage_name, dur in stages.items():
            if dur <= 0:
                continue

            row = {
                "SRR_ID": srr,
                "project_id": project,
                "stage": stage_name,
                "duration_sec": dur,
                "total_reads_M": total_reads / 1e6 if not np.isnan(total_reads) else np.nan,
                "read_length": mean_rl,
                "fastq_mb": fastq_mb,
            }

            if not np.isnan(fastq_mb) and dur > 0:
                row["mb_per_sec"] = fastq_mb / dur
                row["reads_per_sec"] = total_reads / dur if not np.isnan(total_reads) else np.nan
            else:
                row["mb_per_sec"] = np.nan
                row["reads_per_sec"] = np.nan

            # Will be set after category assignment
            row["cores_used"] = 1
            row["mb_per_sec_per_core"] = np.nan

            # Categorize stage
            if stage_name == "fastqc":
                row["category"] = "fastqc"
            elif stage_name.startswith("trim_"):
                row["category"] = "trim"
                row["mode"] = stage_name.replace("trim_", "")
            elif stage_name.startswith("align_"):
                row["category"] = "align"
                row["mode"] = stage_name.replace("align_", "")
            elif stage_name.startswith("count_"):
                row["category"] = "count"
                row["mode"] = stage_name.replace("count_", "")
            elif stage_name.startswith("mode_") and stage_name.endswith("_total"):
                row["category"] = "mode_total"
                row["mode"] = stage_name.replace("mode_", "").replace("_total", "")
            elif stage_name == "pipeline_total":
                row["category"] = "pipeline_total"
            elif stage_name == "sra_conversion":
                row["category"] = "sra_conversion"
            else:
                row["category"] = "other"

            # Set cores and compute per-core throughput
            row["cores_used"] = CORES_PER_STAGE.get(row["category"], 1)
            if not np.isnan(row["mb_per_sec"]):
                row["mb_per_sec_per_core"] = row["mb_per_sec"] / row["cores_used"]

            results.append(row)

    return results


def main():
    parser = argparse.ArgumentParser(description="Calculate MB/s throughput per pipeline stage")
    parser.add_argument("--cores", type=int, default=32, help="Parallel cores")
    parser.add_argument("--project", type=str, default="", help="Process only this project")
    args = parser.parse_args()

    os.makedirs(PLOTS_DIR, exist_ok=True)

    # Load quality data for size estimates (try multiple sources)
    quality_map = {}
    if os.path.exists(QUALITY_FILE):
        with open(QUALITY_FILE) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                srr = row.get("SRR_ID", "")
                if srr:
                    quality_map[srr] = row
        print(f"Loaded quality data for {len(quality_map)} SRRs (from per_srr_quality.tsv)")

    if not quality_map:
        # Fallback 1: extract total reads from featureCounts summaries
        print("Falling back to featureCounts summaries for read counts...")
        quality_map = build_quality_map_from_counts()
        if quality_map:
            print(f"Loaded read counts for {len(quality_map)} SRRs (from featureCounts summaries)")
            print(f"  NOTE: Read length assumed as 150 bp (default). Run 01_extract_fastqc_stats.sh for exact values.")

    added = supplement_quality_from_flattened_bowtie(quality_map)
    if added > 0:
        print(f"Supplemented with {added} SRRs from flattened bowtie stats")

    # Supplement with Bowtie2 alignment stats (legacy aggregated file)
    bowtie_file = os.path.join(OUT_DIR, "bowtie2_alignment_stats.tsv")
    if os.path.exists(bowtie_file):
        added = 0
        with open(bowtie_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                srr = row.get("srr_id", "")
                mode = row.get("mode", "")
                total = safe_float(row.get("total_reads", ""))
                if mode != "untrimmed":
                    continue
                if srr and not np.isnan(total) and srr not in quality_map:
                    quality_map[srr] = {
                        "SRR_ID": srr,
                        "project_id": row.get("project_id", ""),
                        "sequence_depth": str(int(total)),
                        "read_length_mean": "150",  # Assumed default
                    }
                    added += 1
        if added > 0:
            print(f"Supplemented with {added} SRRs from bowtie2_alignment_stats.tsv")

    if not quality_map:
        print(f"WARNING: No quality data found. MB/s calculations will be unavailable.")

    # Discover timing files
    project_dirs = sorted(glob.glob(os.path.join(TIMINGS_BASE, "*/")))
    work_items = []

    for pdir in project_dirs:
        project = os.path.basename(pdir.rstrip("/"))
        if args.project and project != args.project:
            continue
        work_items.append((project, pdir, quality_map))

    print(f"Found {len(work_items)} projects with timing data")

    # Process
    all_results = []
    for item in work_items:
        results = process_project(item)
        all_results.extend(results)

    print(f"Processed {len(all_results)} stage timings")

    if not all_results:
        print("No timing data found.")
        return

    # ============================================================
    # Write Detail Output
    # ============================================================
    detail_file = os.path.join(OUT_DIR, "throughput_detail.tsv")
    cols = ["SRR_ID", "project_id", "category", "stage", "mode",
            "duration_sec", "cores_used", "total_reads_M", "read_length", "fastq_mb",
            "mb_per_sec", "mb_per_sec_per_core", "reads_per_sec"]
    with open(detail_file, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in all_results:
            line = []
            for c in cols:
                v = row.get(c, "")
                if isinstance(v, float):
                    if np.isnan(v):
                        line.append("NA")
                    else:
                        line.append(f"{v:.2f}")
                else:
                    line.append(str(v))
            f.write("\t".join(line) + "\n")
    print(f"Written: {detail_file}")

    # ============================================================
    # Summary by Category and Mode
    # ============================================================
    summary_file = os.path.join(OUT_DIR, "throughput_summary.tsv")
    categories_to_summarize = ["fastqc", "sra_conversion", "trim", "align", "count",
                                "mode_total", "pipeline_total"]

    with open(summary_file, "w") as f:
        f.write("\t".join(["category", "mode", "n", "cores",
                           "dur_mean_sec", "dur_median_sec",
                           "dur_std_sec", "dur_min_sec", "dur_max_sec",
                           "mbps_mean", "mbps_median", "mbps_std",
                           "mbps_per_core_median",
                           "reads_per_sec_mean"]) + "\n")

        for cat in categories_to_summarize:
            # Get unique modes for this category
            modes_in_cat = sorted(set(r.get("mode", "") for r in all_results if r["category"] == cat))
            if not modes_in_cat or modes_in_cat == [""]:
                modes_in_cat = [""]

            for mode in modes_in_cat:
                subset = [r for r in all_results
                          if r["category"] == cat and r.get("mode", "") == mode]
                if not subset:
                    continue

                cores = CORES_PER_STAGE.get(cat, 1)
                durs = np.array([r["duration_sec"] for r in subset])
                mbps = np.array([r["mb_per_sec"] for r in subset
                                 if not np.isnan(r.get("mb_per_sec", np.nan))])
                mbps_pc = np.array([r["mb_per_sec_per_core"] for r in subset
                                    if not np.isnan(r.get("mb_per_sec_per_core", np.nan))])
                rps = np.array([r["reads_per_sec"] for r in subset
                                if not np.isnan(r.get("reads_per_sec", np.nan))])

                f.write("\t".join([
                    cat,
                    mode if mode else "all",
                    str(len(subset)),
                    str(cores),
                    f"{durs.mean():.1f}", f"{np.median(durs):.1f}",
                    f"{durs.std():.1f}", f"{durs.min():.1f}", f"{durs.max():.1f}",
                    f"{mbps.mean():.1f}" if len(mbps) > 0 else "NA",
                    f"{np.median(mbps):.1f}" if len(mbps) > 0 else "NA",
                    f"{mbps.std():.1f}" if len(mbps) > 0 else "NA",
                    f"{np.median(mbps_pc):.2f}" if len(mbps_pc) > 0 else "NA",
                    f"{rps.mean():.0f}" if len(rps) > 0 else "NA",
                ]) + "\n")

    print(f"Written: {summary_file}")

    # ============================================================
    # Print Summary Table
    # ============================================================
    print("\n" + "=" * 100)
    print(f"{'CATEGORY':<20} {'MODE':<15} {'N':>5} {'CORES':>5} {'DUR(s) med':>12} {'MB/s med':>10} {'MB/s/core':>10} {'Reads/s':>12}")
    print("=" * 100)
    for cat in categories_to_summarize:
        modes_in_cat = sorted(set(r.get("mode", "") for r in all_results if r["category"] == cat))
        if not modes_in_cat or modes_in_cat == [""]:
            modes_in_cat = [""]
        for mode in modes_in_cat:
            subset = [r for r in all_results
                      if r["category"] == cat and r.get("mode", "") == mode]
            if not subset:
                continue
            cores = CORES_PER_STAGE.get(cat, 1)
            durs = np.array([r["duration_sec"] for r in subset])
            mbps = [r["mb_per_sec"] for r in subset if not np.isnan(r.get("mb_per_sec", np.nan))]
            mbps_pc = [r["mb_per_sec_per_core"] for r in subset if not np.isnan(r.get("mb_per_sec_per_core", np.nan))]
            rps = [r["reads_per_sec"] for r in subset if not np.isnan(r.get("reads_per_sec", np.nan))]
            mbps_str = f"{np.median(mbps):>10.1f}" if mbps else f"{'NA':>10}"
            mbpc_str = f"{np.median(mbps_pc):>10.2f}" if mbps_pc else f"{'NA':>10}"
            rps_str = f"{np.mean(rps):>11.0f}" if rps else f"{'NA':>11}"
            print(f"{cat:<20} {mode or 'all':<15} {len(subset):>5} {cores:>5} "
                  f"{np.median(durs):>12.1f} {mbps_str} {mbpc_str}  {rps_str}")
    print("=" * 100)

    # ============================================================
    # Plots
    # ============================================================
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Skipping plots (matplotlib not available)")
        return

    # --- Plot 1: Core-seconds boxplots by stage ---
    stage_order = ["fastqc", "sra_conversion"]
    for mode in ALIGN_MODES:
        for prefix in ["trim", "align", "count"]:
            stage_order.append(f"{prefix}_{mode}")

    fig, ax = plt.subplots(figsize=(20, 8))
    data_for_plot = []
    labels = []
    for stage in stage_order:
        stage_results = [r for r in all_results if r["stage"] == stage]
        if stage_results:
            core_secs = [r["duration_sec"] * r["cores_used"] for r in stage_results]
            data_for_plot.append(core_secs)
            labels.append(stage.replace("_", "\n"))

    if data_for_plot:
        bp = ax.boxplot(data_for_plot, labels=labels, patch_artist=True, notch=True,
                        widths=0.6, flierprops=dict(markersize=3, alpha=0.4))
        # Color by category
        cat_colors = {"fastqc": "#42A5F5", "sra": "#78909C",
                      "trim": "#66BB6A", "align": "#FFA726", "count": "#AB47BC"}
        for i, label_text in enumerate(labels):
            flat = label_text.replace("\n", "_")
            if flat.startswith("trim"):
                bp["boxes"][i].set_facecolor(cat_colors["trim"])
            elif flat.startswith("align"):
                bp["boxes"][i].set_facecolor(cat_colors["align"])
            elif flat.startswith("count"):
                bp["boxes"][i].set_facecolor(cat_colors["count"])
            elif flat.startswith("fastqc"):
                bp["boxes"][i].set_facecolor(cat_colors["fastqc"])
            else:
                bp["boxes"][i].set_facecolor(cat_colors["sra"])
            bp["boxes"][i].set_alpha(0.7)

        ax.set_ylabel("Core-seconds (duration × cores)", fontsize=12)
        ax.set_title("Computational Cost per Pipeline Stage (All SRRs)", fontsize=14)
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.grid(axis="y", alpha=0.3)
        ax.set_yscale("log")

        # Legend
        from matplotlib.patches import Patch
        legend_elems = [Patch(facecolor=c, alpha=0.7, label=f"{l} ({CORES_PER_STAGE.get(l.lower().rstrip('ming').rstrip('ment'), '?')}c)")
                        for l, c in [("FastQC", "#42A5F5"), ("Trimming", "#66BB6A"),
                                     ("Alignment", "#FFA726"), ("Counting", "#AB47BC")]]
        legend_elems = [
            Patch(facecolor="#42A5F5", alpha=0.7, label="FastQC (8 cores)"),
            Patch(facecolor="#66BB6A", alpha=0.7, label="Trimming (36 cores)"),
            Patch(facecolor="#FFA726", alpha=0.7, label="Alignment (36 cores)"),
            Patch(facecolor="#AB47BC", alpha=0.7, label="Counting (36 cores)"),
        ]
        ax.legend(handles=legend_elems, loc="upper right")

    fig.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, "stage_core_seconds.png")
    fig.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Written: {plot_file}")

    # --- Plot 2: MB/s per core throughput by category ---
    fig, axes = plt.subplots(1, 4, figsize=(18, 6), sharey=False)
    cat_labels = ["fastqc", "trim", "align", "count"]
    cat_colors_list = ["#42A5F5", "#66BB6A", "#FFA726", "#AB47BC"]

    for ax_i, (cat, color) in enumerate(zip(cat_labels, cat_colors_list)):
        ax = axes[ax_i]
        cores = CORES_PER_STAGE.get(cat, 1)
        modes = sorted(set(r.get("mode", "all") for r in all_results if r["category"] == cat))
        if not modes or modes == [""]:
            modes = ["all"]
        data = []
        mode_labels = []
        for mode in modes:
            vals = [r["mb_per_sec_per_core"] for r in all_results
                    if r["category"] == cat and r.get("mode", "all") == mode
                    and not np.isnan(r.get("mb_per_sec_per_core", np.nan))]
            if vals:
                data.append(vals)
                mode_labels.append(mode if mode else "all")

        if data:
            bp = ax.boxplot(data, labels=mode_labels, patch_artist=True, notch=True,
                            widths=0.6, flierprops=dict(markersize=3, alpha=0.4))
            for patch in bp["boxes"]:
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
        ax.set_title(f"{cat.upper()} ({cores} cores)", fontsize=12)
        ax.set_ylabel("MB/s per core" if ax_i == 0 else "", fontsize=11)
        ax.tick_params(axis="x", rotation=45, labelsize=9)
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Throughput per Core (MB/s/core) by Pipeline Stage", fontsize=14, y=1.02)
    fig.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, "throughput_per_core.png")
    fig.savefig(plot_file, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Written: {plot_file}")

    # --- Plot 3: Core-seconds vs Read Count scatter ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    scatter_cats = [("fastqc", "#42A5F5"), ("trim", "#66BB6A"),
                    ("align", "#FFA726"), ("count", "#AB47BC")]
    for ax_i, (cat, color) in enumerate(scatter_cats):
        ax = axes[ax_i // 2][ax_i % 2]
        cores = CORES_PER_STAGE.get(cat, 1)
        subset = [r for r in all_results if r["category"] == cat
                  and not np.isnan(r.get("total_reads_M", np.nan))]
        if subset:
            x = [r["total_reads_M"] for r in subset]
            y = [r["duration_sec"] * r["cores_used"] for r in subset]
            ax.scatter(x, y, c=color, alpha=0.4, s=15, edgecolors="none")
            # Fit trend line
            if len(x) > 5:
                z = np.polyfit(x, y, 1)
                x_fit = np.linspace(min(x), max(x), 100)
                ax.plot(x_fit, np.polyval(z, x_fit), "--", color="black", alpha=0.5,
                        label=f"slope={z[0]:.1f} core-s/M reads")
                ax.legend(fontsize=9)
        ax.set_xlabel("Total Reads (millions)", fontsize=10)
        ax.set_ylabel("Core-seconds", fontsize=10)
        ax.set_title(f"{cat.upper()} ({cores} cores)", fontsize=12)
        ax.grid(alpha=0.3)

    fig.suptitle("Computational Cost vs Sequencing Depth", fontsize=14, y=1.02)
    fig.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, "core_seconds_vs_reads.png")
    fig.savefig(plot_file, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Written: {plot_file}")

    print("\nThroughput analysis complete.")


if __name__ == "__main__":
    main()
