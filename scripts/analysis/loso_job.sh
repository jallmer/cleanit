#!/bin/bash
#SBATCH --job-name=loso_gsea
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/logs/loso_%j.out
#SBATCH --error=/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/logs/loso_%j.err

PROJECT_ID=$1
# Default to SLURM cores if 2nd arg is not provided
CORES=${2:-${SLURM_CPUS_PER_TASK:-16}}

echo "Starting LOSO for project: $PROJECT_ID"

# Prevent numpy/OpenBLAS from spawning 16 threads per python worker (16 workers * 16 threads = 256 threads)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

source /pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/bin/activate omiks_pipeline
python3 /pc2/users/o/omiks001/scripts/analysis/06_run_de_gsea.py --project "$PROJECT_ID" --cores "$CORES"
echo "Finished LOSO for project: $PROJECT_ID"
