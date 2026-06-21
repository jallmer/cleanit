#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────
#  install_sexcall.sh – Set up conda environment for sexcall
# ─────────────────────────────────────────────────────────────────────
#  Creates the env_pysam environment with pysam (and samtools).
#  Safe to re-run – skips if the environment already exists.
# ─────────────────────────────────────────────────────────────────────

set -euo pipefail

export CONDA_SOLVER=classic
export CONDA_NO_PLUGINS=true

# ── Channel hygiene ─────────────────────────────────────────────────
echo ">> Setting channel order (bioconda → conda-forge → defaults) ..."
for c in bioconda conda-forge defaults; do
    conda config --remove channels "$c" >/dev/null 2>&1 || true
done
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels defaults
conda config --set channel_priority strict

# ── Create env_pysam ────────────────────────────────────────────────
ENV_NAME="env_pysam"
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo ">> $ENV_NAME already exists – skipping."
else
    echo ">> Creating $ENV_NAME ..."
    conda create -y -n "$ENV_NAME" bioconda::pysam bioconda::samtools \
          -c bioconda -c conda-forge -c defaults
fi

# ── Make scripts executable ─────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
chmod +x "$SCRIPT_DIR/sexcall.sh" "$SCRIPT_DIR/run_sexcall_batch.sh" 2>/dev/null || true

cat <<EOF

────────────────────────────────────────────────────────────────
  ✔  sexcall is ready.

  Quick start:
    ./sexcall.sh  <sorted.bam>
    ./sexcall.sh  <sorted.bam>  --xist XIST_exons.hg38.bed

  Batch mode:
    ./run_sexcall_batch.sh  --xist XIST_exons.hg38.bed

  See ./sexcall.sh --help for all options.
────────────────────────────────────────────────────────────────
EOF
