#!/bin/bash
# v2_pure_submitter.sh
# Usage: bash v2_pure_submitter.sh [NUM_JOBS]

DB_FILE="$HOME/srr_queue.db"
NUM_JOBS=${1:-50}

echo "Initializing Stage 1 Sequence Mapping ('No Bullshit Mode')."
echo "Extracting exactly up to $NUM_JOBS queued sequences..."

JOBS=$(sqlite3 "$DB_FILE" "SELECT project_id, srr_id, local_fastq_path FROM srr_queue WHERE status='todo' AND run_status='FETCHED' AND local_fastq_path IS NOT NULL LIMIT $NUM_JOBS;")

if [ -z "$JOBS" ]; then
    echo "Database formally empty of FETCHED bounds. Zero completely pending natively found. Run v2_data_fetcher.sh first!"
    exit 0
fi

for job in $JOBS; do
    PRJ=$(echo "$job" | cut -d '|' -f 1)
    SRR=$(echo "$job" | cut -d '|' -f 2)
    FASTQ_P=$(echo "$job" | cut -d '|' -f 3)
    
    echo "Queueing monolith pipeline structure strictly for $SRR ($PRJ)..."
    
    # Fire identical master script algorithm into SLURM queue physically!
    # Inherently coupling 12 Cores tightly to 48GB structurally to match the user's Billing Profile dynamically!
    sbatch /pc2/users/o/omiks001/scripts/v2_master_node.sh "$PRJ" "$SRR" "PAIRED" "$FASTQ_P"
    
    # Save the queue bounds directly
    sqlite3 "$DB_FILE" "UPDATE srr_queue SET status='submitted', run_status='SLURM QUEUED' WHERE srr_id='$SRR';"
    
    # Protect SLURM control node securely against instantaneous queue blocking
    sleep 0.5
done

echo "Submission sequence computationally terminated securely!"
