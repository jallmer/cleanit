#!/bin/bash
# 03_extract_trimmomatic_stats.sh — Aggregate trimming and bowtie stats for analysis.
# Usage: bash scripts/analysis/03_extract_trimmomatic_stats.sh
#
# Preferred inputs:
#   /scratch/hpc-prf-omiks/ja/flattened_trimmomatic_stats/<PROJECT>/<SRR>_trimmomatic_stats.tsv
#   /scratch/hpc-prf-omiks/ja/flattened_bowtie2_stats/<PROJECT>/<SRR>_bowtie2_stats.tsv
#
# Legacy fallback:
#   ~/logs/v2_pipeline_*.out
#   ~/logs/v2_pipeline_*.err
#
# Writes to:
#   /scratch/hpc-prf-omiks/ja/analysis/trimmomatic_detail.tsv
#   /scratch/hpc-prf-omiks/ja/analysis/trimmomatic_summary_by_mode.tsv
#   /scratch/hpc-prf-omiks/ja/analysis/trimmomatic_summary_by_project.tsv
#   /scratch/hpc-prf-omiks/ja/analysis/bowtie2_alignment_stats.tsv

set -euo pipefail

OUT_DIR="/scratch/hpc-prf-omiks/ja/analysis"
DETAIL_FILE="$OUT_DIR/trimmomatic_detail.tsv"
SUMMARY_MODE="$OUT_DIR/trimmomatic_summary_by_mode.tsv"
SUMMARY_PROJ="$OUT_DIR/trimmomatic_summary_by_project.tsv"
BOWTIE_FILE="$OUT_DIR/bowtie2_alignment_stats.tsv"
TRIM_BASE="/scratch/hpc-prf-omiks/ja/flattened_trimmomatic_stats"
BOWTIE_BASE="/scratch/hpc-prf-omiks/ja/flattened_bowtie2_stats"

mkdir -p "$OUT_DIR"

echo "============================================"
echo "Flattened Stats Aggregation"
echo "============================================"
echo "Trimmomatic source: $TRIM_BASE"
echo "Bowtie2 source:     $BOWTIE_BASE"
echo "Output dir:         $OUT_DIR"
echo ""

python3 - "$TRIM_BASE" "$BOWTIE_BASE" "$DETAIL_FILE" "$SUMMARY_MODE" "$SUMMARY_PROJ" "$BOWTIE_FILE" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

trim_base = Path(sys.argv[1])
bowtie_base = Path(sys.argv[2])
detail_file = Path(sys.argv[3])
summary_mode = Path(sys.argv[4])
summary_proj = Path(sys.argv[5])
bowtie_file = Path(sys.argv[6])

trim_header = ["project_id", "SRR_ID", "mode", "layout", "input_reads", "surviving", "surviving_pct", "dropped", "dropped_pct"]
bowtie_header = ["job_id", "srr_id", "project_id", "mode", "layout", "total_reads", "paired_reads", "unpaired_reads",
                 "aligned_exactly_1", "aligned_gt1", "aligned_0", "concordant_exactly_1", "concordant_gt1",
                 "discordant_1", "pairs_0_concordant", "overall_alignment_rate"]

trim_rows = []
if trim_base.exists():
    for path in sorted(trim_base.glob("*/*.tsv")):
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                trim_rows.append({
                    "project_id": row.get("project_id", ""),
                    "SRR_ID": row.get("srr_id", row.get("SRR_ID", "")),
                    "mode": row.get("mode", ""),
                    "layout": row.get("layout", ""),
                    "input_reads": row.get("input_reads", ""),
                    "surviving": row.get("surviving", ""),
                    "surviving_pct": row.get("surviving_pct", ""),
                    "dropped": row.get("dropped", ""),
                    "dropped_pct": row.get("dropped_pct", ""),
                })

with detail_file.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=trim_header, delimiter="\t")
    writer.writeheader()
    writer.writerows(trim_rows)

def to_int(v):
    try:
        return int(v)
    except Exception:
        return 0

def to_float(v):
    try:
        return float(str(v).rstrip("%"))
    except Exception:
        return 0.0

by_mode = defaultdict(lambda: {"samples": 0, "input": 0, "surv": 0, "pct_sum": 0.0})
by_proj_mode = defaultdict(lambda: {"samples": 0, "pct_sum": 0.0})
for row in trim_rows:
    mode = row["mode"]
    key = (row["project_id"], mode)
    by_mode[mode]["samples"] += 1
    by_mode[mode]["input"] += to_int(row["input_reads"])
    by_mode[mode]["surv"] += to_int(row["surviving"])
    by_mode[mode]["pct_sum"] += to_float(row["surviving_pct"])
    by_proj_mode[key]["samples"] += 1
    by_proj_mode[key]["pct_sum"] += to_float(row["surviving_pct"])

with summary_mode.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["mode", "samples", "total_input", "total_surviving", "surviving_pct", "mean_surviving_pct", "total_dropped", "dropped_pct"])
    for mode in sorted(by_mode):
        rec = by_mode[mode]
        dropped = rec["input"] - rec["surv"]
        surv_pct = (rec["surv"] / rec["input"] * 100) if rec["input"] else 0.0
        drop_pct = (dropped / rec["input"] * 100) if rec["input"] else 0.0
        mean_pct = rec["pct_sum"] / rec["samples"] if rec["samples"] else 0.0
        writer.writerow([mode, rec["samples"], rec["input"], rec["surv"], f"{surv_pct:.2f}%", f"{mean_pct:.2f}%", dropped, f"{drop_pct:.2f}%"])

with summary_proj.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["project_id", "mode", "samples", "mean_surviving_pct"])
    for (project_id, mode) in sorted(by_proj_mode):
        rec = by_proj_mode[(project_id, mode)]
        mean_pct = rec["pct_sum"] / rec["samples"] if rec["samples"] else 0.0
        writer.writerow([project_id, mode, rec["samples"], f"{mean_pct:.2f}%"])

bowtie_rows = []
if bowtie_base.exists():
    for path in sorted(bowtie_base.glob("*/*.tsv")):
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                bowtie_rows.append({
                    "job_id": "NA",
                    "srr_id": row.get("srr_id", ""),
                    "project_id": row.get("project_id", ""),
                    "mode": row.get("mode", ""),
                    "layout": row.get("layout", ""),
                    "total_reads": row.get("total_reads", ""),
                    "paired_reads": row.get("paired_reads", ""),
                    "unpaired_reads": row.get("unpaired_reads", ""),
                    "aligned_exactly_1": row.get("aligned_exactly_1", ""),
                    "aligned_gt1": row.get("aligned_gt1", ""),
                    "aligned_0": row.get("aligned_0", ""),
                    "concordant_exactly_1": row.get("concordant_exactly_1", ""),
                    "concordant_gt1": row.get("concordant_gt1", ""),
                    "discordant_1": row.get("discordant_1", ""),
                    "pairs_0_concordant": row.get("pairs_0_concordant", ""),
                    "overall_alignment_rate": row.get("overall_alignment_rate", ""),
                })

with bowtie_file.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=bowtie_header, delimiter="\t")
    writer.writeheader()
    writer.writerows(bowtie_rows)

print(f"trimmomatic_rows\t{len(trim_rows)}")
print(f"bowtie_rows\t{len(bowtie_rows)}")
PY

echo ""
echo "Written:"
echo "  $DETAIL_FILE"
echo "  $SUMMARY_MODE"
echo "  $SUMMARY_PROJ"
echo "  $BOWTIE_FILE"
