#!/bin/bash
# 05_run_all.sh — Master runner for the analysis suite
# Usage: bash scripts/analysis/05_run_all.sh [PROJECT_FILTER] [CORES]
#
# Runs all analysis steps in sequence:
#   00: Environment check
#   01: FastQC quality extraction
#   02: Count matrix evaluation (JSD + Pearson)
#   03: Trimmomatic stats extraction
#   04: Summary aggregation + plots
#
# All output goes to /scratch/hpc-prf-omiks/ja/analysis/

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="${1:-}"
CORES="${2:-32}"

echo "============================================================"
echo "   RNA-Seq Analysis Suite — Master Runner"
echo "============================================================"
echo "Directory: $DIR"
echo "Cores:     $CORES"
[ -n "$PROJECT" ] && echo "Project:   $PROJECT"
echo "Start:     $(date)"
echo "============================================================"
echo ""

# Activate conda if needed
if ! command -v python3 &>/dev/null; then
    echo "Activating conda..."
    source /scratch/hpc-prf-omiks/ja/miniconda3/etc/profile.d/conda.sh
    conda activate base
fi

TOTAL_START=$SECONDS

# Step 0: Environment check
echo ">>> STEP 0: Environment Check"
python3 -c "import numpy; import scipy" 2>/dev/null || {
    echo "Installing dependencies..."
    bash "$DIR/00_setup_env.sh"
}
echo "  Python: $(python3 --version)"
echo "  numpy:  $(python3 -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo 'MISSING')"
echo "  scipy:  $(python3 -c 'import scipy; print(scipy.__version__)' 2>/dev/null || echo 'MISSING')"
echo ""

# Step 1: FastQC quality extraction
echo ">>> STEP 1: FastQC Quality Extraction"
S1=$SECONDS
bash "$DIR/01_extract_fastqc_stats.sh" "$PROJECT"
echo "  Time: $(( SECONDS - S1 ))s"
echo ""

# Step 2: Count matrix evaluation
echo ">>> STEP 2: Count Matrix Evaluation (JSD + Pearson)"
S2=$SECONDS
if [ -n "$PROJECT" ]; then
    python3 "$DIR/02_evaluate_count_matrices.py" --cores "$CORES" --project "$PROJECT"
else
    python3 "$DIR/02_evaluate_count_matrices.py" --cores "$CORES"
fi
echo "  Time: $(( SECONDS - S2 ))s"
echo ""

# Step 3: Trimmomatic stats (only if logs exist)
echo ">>> STEP 3: Flattened Stats Aggregation"
S3=$SECONDS
bash "$DIR/03_extract_trimmomatic_stats.sh"
echo "  Time: $(( SECONDS - S3 ))s"
echo ""


# Step 4: Summary aggregation + plots
echo ">>> STEP 4: Summary Aggregation + Plots"
S4=$SECONDS
python3 "$DIR/04_summarize_results.py"
echo "  Time: $(( SECONDS - S4 ))s"
echo ""

# Step 5: Throughput
echo ">>> STEP 5: Throughput Analysis"
S5=$SECONDS
if [ -n "$PROJECT" ]; then
    python3 "$DIR/11_compute_throughput.py" --cores "$CORES" --project "$PROJECT"
else
    python3 "$DIR/11_compute_throughput.py" --cores "$CORES"
fi
echo "  Time: $(( SECONDS - S5 ))s"
echo ""

TIER1_DUR=$(( SECONDS - TOTAL_START ))
echo "============================================================"
echo "   TIER 1 COMPLETE  (${TIER1_DUR}s)"
echo "============================================================"
echo ""

# ============================================================
# TIER 2: DE/GSEA Methodology (requires condition metadata)
# ============================================================
METADATA_DIR="$DIR/metadata"
HAS_METADATA=false

if [ -d "$METADATA_DIR" ] && ls "$METADATA_DIR"/*.csv &>/dev/null; then
    META_COUNT=$(ls "$METADATA_DIR"/*.csv 2>/dev/null | grep -cv "_attributes" || true)
    if [ "$META_COUNT" -gt 0 ]; then
        HAS_METADATA=true
        echo "Found $META_COUNT project metadata files in $METADATA_DIR"
    fi
fi

if [ "$HAS_METADATA" = true ]; then
    # Step 6: DE + GSEA concordance
    echo ""
    echo ">>> STEP 6: Leave-One-Out DE + GSEA Concordance"
    S6=$SECONDS
    if [ -n "$PROJECT" ]; then
        python3 "$DIR/06_run_de_gsea.py" --cores "$CORES" --project "$PROJECT" --metadata-dir "$METADATA_DIR"
    else
        python3 "$DIR/06_run_de_gsea.py" --cores "$CORES" --metadata-dir "$METADATA_DIR"
    fi
    echo "  Time: $(( SECONDS - S6 ))s"
    echo ""

    # Step 7: Classification
    echo ">>> STEP 7: Trimming Classification"
    S7=$SECONDS
    python3 "$DIR/07_classify_trimming.py"
    echo "  Time: $(( SECONDS - S7 ))s"
    echo ""

    # Step 8: Build merged feature table
    echo ">>> STEP 8: Build Sample Feature Table"
    S8=$SECONDS
    python3 "$DIR/13_build_sample_feature_table.py"
    echo "  Time: $(( SECONDS - S8 ))s"
    echo ""

    # Step 9: QC Model
    echo ">>> STEP 9: QC Predictor Model"
    S9=$SECONDS
    python3 "$DIR/08_fit_qc_model.py"
    echo "  Time: $(( SECONDS - S9 ))s"
    echo ""

    # Step 10: Reports
    echo ">>> STEP 10: Generate Reports"
    S10=$SECONDS
    python3 "$DIR/09_generate_reports.py"
    echo "  Time: $(( SECONDS - S10 ))s"
    echo ""
else
    echo "No condition metadata found in $METADATA_DIR"
    echo "  To enable Tier 2 analysis (DE/GSEA), place Run,condition CSVs there."
    echo "  Use: python3 $DIR/10_fetch_metadata.py --project PRJNA..."
    echo ""
fi

TOTAL_DUR=$(( SECONDS - TOTAL_START ))
echo "============================================================"
echo "   ALL DONE  (Total: ${TOTAL_DUR}s)"
echo "============================================================"
echo ""
echo "Output directory: /scratch/hpc-prf-omiks/ja/analysis/"
echo ""
echo "Key files (Tier 1):"
echo "  per_srr_quality.tsv        — FastQC quality metrics"
echo "  per_srr_eval.tsv           — JSD + Pearson per SRR"
echo "  global_summary.tsv         — Merged quality + eval"
echo "  per_project_summary.tsv    — Per-project aggregated stats"
echo "  per_project/<PROJ>_eval.tsv — Per-project eval detail"
echo "  trimmomatic_detail.tsv     — Trimming stats"
echo "  bowtie2_alignment_stats.tsv — Bowtie2 stats"
echo "  throughput_detail.tsv      — Stage throughput detail"
echo "  throughput_summary.tsv     — Stage throughput summary"
echo "  plots/                     — Boxplots and heatmaps"
echo "  multiqc/                   — MultiQC reports (if available)"
if [ "$HAS_METADATA" = true ]; then
    echo ""
    echo "Key files (Tier 2 — DE/GSEA):"
    echo "  concordance/               — Per-project concordance tables"
    echo "  trimming_classification.tsv — Helpful/neutral/harmful per SRR"
    echo "  trimming_benefit.tsv       — Benefit scores per SRR"
    echo "  trimming_penalties.tsv     — Read/gene/assigned-fraction losses per SRR-method"
    echo "  sample_feature_table.tsv   — Merged QC + penalty + benefit feature table"
    echo "  qc_model_results.tsv       — Mixed-effects model output"
    echo "  qc_model_optimal_method.tsv — Ordinal model for t*"
    echo "  classification_summary.tsv — Summary by project"
fi
