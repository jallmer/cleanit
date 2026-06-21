#!/bin/bash
#SBATCH --job-name=v2_supervisor
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --partition=normal
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G

echo "=========================================================="
echo "Starting Native V2 Pipeline Supervisor Daemon via SLURM"
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo "=========================================================="

# Run the supervisor daemon structurally natively inside the cluster boundary!
bash /pc2/users/o/omiks001/scripts/v2_supervisor.sh
