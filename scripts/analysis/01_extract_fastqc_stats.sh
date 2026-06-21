#!/bin/bash
# 01_extract_fastqc_stats.sh — Extract quality metrics from FastQC output directories
# Usage: bash scripts/analysis/01_extract_fastqc_stats.sh [PROJECT_FILTER]
# Reads from: /scratch/hpc-prf-omiks/ja/flattened_fastqc_raw/
# Writes to:  /scratch/hpc-prf-omiks/ja/analysis/per_srr_quality.tsv

set -euo pipefail

FASTQC_BASE="/scratch/hpc-prf-omiks/ja/flattened_fastqc_raw"
OUT_DIR="/scratch/hpc-prf-omiks/ja/analysis"
OUT_FILE="$OUT_DIR/per_srr_quality.tsv"
TMP_FILE="$OUT_DIR/per_srr_quality.raw.tsv"
PROJECT_FILTER="${1:-}"
CORES="${2:-32}"

echo "============================================"
echo "FastQC Quality Extraction"
echo "============================================"
echo "Source: $FASTQC_BASE"
echo "Output: $OUT_FILE"
[ -n "$PROJECT_FILTER" ] && echo "Filter: $PROJECT_FILTER"
echo ""

python3 "$HOME/scripts/analysis/01_extract_fastqc_stats.py" "$PROJECT_FILTER"

TOTAL=$(tail -n +2 "$OUT_FILE" | wc -l)
UNIQUE_SRRS=$(tail -n +2 "$OUT_FILE" | cut -f2 | sort -u | wc -l)
echo ""
echo "Done. Total FastQC entries: $TOTAL, Unique SRRs: $UNIQUE_SRRS"
echo "Output: $OUT_FILE"

# Run MultiQC per project if available
if command -v multiqc &>/dev/null; then
    echo ""
    echo "Running MultiQC per project..."
    for project_dir in "$FASTQC_BASE"/*/; do
        project=$(basename "$project_dir")
        [ -n "$PROJECT_FILTER" ] && [ "$project" != "$PROJECT_FILTER" ] && continue
        mqc_out="$OUT_DIR/multiqc/${project}"
        if [ ! -f "$mqc_out/multiqc_report.html" ]; then
            echo "  MultiQC: $project"
            multiqc -q -o "$mqc_out" "$project_dir" 2>/dev/null || echo "    WARN: MultiQC failed for $project"
        fi
    done
else
    echo "MultiQC not found — skipping aggregated reports."
fi
