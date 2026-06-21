#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/deseq_%x_%j.out
#SBATCH --error=logs/deseq_%x_%j.err

# Run whole-project DESeq2 + GSEA for a SINGLE project.
# One core per method (U, A, P5, P10, P20, P35), 4GB per core.
#
# Usage:
#   sbatch --job-name=PRJNA1105191 scripts/submit_deseq2_gsea.sh PRJNA1105191

ANALYSIS="/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis"
PYTHON="/scratch/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/python"

PROJECT="$1"
if [ -z "$PROJECT" ]; then
    echo "ERROR: No project ID provided"
    exit 1
fi

cd "$ANALYSIS"

echo "Project: $PROJECT"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Mem per CPU: 4G"
echo "Started: $(date)"
echo ""

$PYTHON scripts/run_deseq2_gsea.py "$PROJECT"

echo ""
echo "Finished: $(date)"
