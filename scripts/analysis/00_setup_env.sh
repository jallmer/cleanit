#!/bin/bash
# 00_setup_env.sh — Install Python packages needed for the analysis suite
# Run once: bash scripts/analysis/00_setup_env.sh

set -euo pipefail

echo "============================================"
echo "Analysis Environment Setup"
echo "============================================"

# Output directory
OUT_DIR="/scratch/hpc-prf-omiks/ja/analysis"
mkdir -p "$OUT_DIR"/{per_project,plots,multiqc}

# Check Python
if ! command -v python &>/dev/null; then
    echo "ERROR: python not found. Activate conda first:"
    echo "  source /pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/etc/profile.d/conda.sh"
    echo "  conda activate base"
    exit 1
fi

echo "Python: $(which python) — $(python --version)"

# Install minimal deps via pip (into user or conda env)
echo ""
echo "Installing Python dependencies..."
pip install --quiet numpy scipy matplotlib multiqc pydeseq2 gseapy statsmodels anndata 2>/dev/null || \
    conda install -y -c conda-forge -c bioconda numpy scipy matplotlib multiqc pydeseq2 gseapy statsmodels anndata 2>/dev/null || \
    echo "WARNING: Could not install some packages. Verify manually."

# Verify
echo ""
echo "Checking dependencies..."
python -c "import numpy; print(f'  numpy {numpy.__version__}')"
python -c "import scipy; print(f'  scipy {scipy.__version__}')"
python -c "import matplotlib; print(f'  matplotlib {matplotlib.__version__}')"
command -v multiqc && echo "  multiqc $(multiqc --version 2>&1 | head -1)" || echo "  multiqc: NOT FOUND (optional)"

echo ""
echo "Output directory: $OUT_DIR"
echo "Setup complete."
