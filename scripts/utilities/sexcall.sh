#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────
#  sexcall.sh – Wrapper for sexcall.py
# ─────────────────────────────────────────────────────────────────────
#  Usage:
#    sexcall.sh <sorted.bam> [--xist XIST.bed] [--mapq 30] [--y-threshold 10]
#                            [--output out.tsv] [--append] [--sample NAME]
#
#  Activates the pysam conda env and invokes the Python caller.
# ─────────────────────────────────────────────────────────────────────

set -euo pipefail

# Location of the Python script (same directory as this wrapper)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEXCALL_PY="$SCRIPT_DIR/sexcall.py"

if [[ ! -f "$SEXCALL_PY" ]]; then
    echo "[ERROR] sexcall.py not found at $SEXCALL_PY" >&2
    exit 1
fi

if [[ $# -lt 1 || "$1" == "-h" || "$1" == "--help" ]]; then
    python3 "$SEXCALL_PY" --help
    exit 0
fi

# ── Activate conda environment ──────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_pysam

# ── Run the caller ──────────────────────────────────────────────────
python3 "$SEXCALL_PY" "$@"
