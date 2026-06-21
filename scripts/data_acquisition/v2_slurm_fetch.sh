#!/bin/bash
#SBATCH --job-name=v2_fetch
#SBATCH --output=/pc2/users/o/omiks001/scripts/fetch_logs/%x_%j.log
#SBATCH --error=/pc2/users/o/omiks001/scripts/fetch_logs/%x_%j.err
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Arguments passed via: sbatch v2_slurm_fetch.sh <PROJECT> <SRR> <ENDPOINT>
PROJECT=$1
SRR=$2
ENDPOINT=$3

SCRIPT_DIR="/pc2/users/o/omiks001/scripts"
CONDA_SH="/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/etc/profile.d/conda.sh"
CONDA_BASE_BIN="/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/bin"
PIPELINE_ENV_BIN="/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin"

if [ -f "$CONDA_SH" ]; then
    . "$CONDA_SH"
fi
export PATH="$PIPELINE_ENV_BIN:$CONDA_BASE_BIN:$PATH"
if command -v conda >/dev/null 2>&1; then
    conda activate omiks_pipeline >/dev/null 2>&1 || true
fi

for tool in prefetch ascp curl wget; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: required fetch tool '$tool' is not available on PATH"
        exit 1
    fi
done

echo "Fetch job started: $SRR via $ENDPOINT on $(hostname) at $(date)"

bash "$SCRIPT_DIR/v2_fetch_worker.sh" "$PROJECT" "$SRR" "$ENDPOINT"
EXIT_CODE=$?

echo "Fetch job finished: $SRR exit=$EXIT_CODE at $(date)"
exit $EXIT_CODE
