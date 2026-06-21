#!/bin/bash
#SBATCH -J v2_macro_array
#SBATCH -o macro_%j.out
#SBATCH -e macro_%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=36  # Divisible by 6 perfectly
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH -p normal

# The Macro Runner mathematically aggregates multiple SQLite jobs into a single persistent SLURM allocation!
# This fundamentally bypasses the 1-job-per-sequence scaling limitation structurally.

NUM_JOBS=${1:-20}
DB_FILE="$HOME/srr_queue.db"
export SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

echo "Macro Aggregator Starting... Pulling exactly $NUM_JOBS sequences computationally."

# Extract structurally available sequences
JOBS=$(sqlite3 "$DB_FILE" "SELECT project_id, srr_id FROM srr_queue WHERE status='todo' LIMIT $NUM_JOBS;")

if [ -z "$JOBS" ]; then
    echo "No todo sequences found in the registry."
    exit 0
fi

for job in $JOBS; do
    PRJ=$(echo "$job" | cut -d '|' -f 1)
    SRR=$(echo "$job" | cut -d '|' -f 2)
    
    echo "=========================================================="
    echo "MACRO NODE: Initiating mathematical wrapper for $SRR"
    echo "=========================================================="
    
    # Secure ownership in the global database natively
    sqlite3 "$DB_FILE" "UPDATE srr_queue SET status='submitted', run_status='INITIALIZED', fetch_lane='MACRO', slurm_job_id='$SLURM_JOB_ID' WHERE srr_id='$SRR';"
    
    # Execute exactly one node layer explicitly natively through the universal sequence mapper (NCBI fallbacks internally!)
    bash "/pc2/users/o/omiks001/scripts/v2_master_node.sh" "$PRJ" "$SRR" ncbi
    
    EXIT_CODE=$?
    if [ $EXIT_CODE -ne 0 ]; then
        echo "MACRO NODE: Sequence $SRR threw a fatal pipeline exit. Continuing to next sequence..." >&2
    fi
    
    # Sleep to allow memory buffers to physically clear gracefully
    sleep 10
done

echo "Macro Array Execution Completed Successfully for $NUM_JOBS Sequences!"
