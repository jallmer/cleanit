#!/bin/bash
# v2_submitter.sh
# Usage: bash v2_submitter.sh <LIMIT> [FETCH_MODE]
# Reads precisely unsubmitted SRR records from the SQLite tracker, blocks them off to 'submitted', and natively SLURMs the monolithic framework.

LIMIT=${1:-10}
DB_FILE="$HOME/srr_queue.db"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python3 "$DIR/v2_reconcile_queue.py" --db "$DB_FILE" --apply >/dev/null 2>&1 || true

echo "Fetching up to $LIMIT 'FETCHED' computation matrices directly off the SQL table -> $DB_FILE..."

# Note: sqlite3 natively parses identically strings structurally!
JOBS=$(sqlite3 "$DB_FILE" "SELECT project_id, srr_id, local_fastq_path
FROM srr_queue
WHERE status = 'todo'
  AND run_status = 'FETCHED'
  AND local_fastq_path IS NOT NULL
  AND local_fastq_path != ''
ORDER BY failure_count ASC LIMIT $LIMIT;")

if [ -z "$JOBS" ]; then
    echo "No 'FETCHED' arrays formally scraped globally. Trigger v2_data_fetcher.sh locally!"
    exit 0
fi

while IFS='|' read -r PROJECT_ID SRR_ID FASTQ_PATH; do
    [ -z "$PROJECT_ID" ] && continue
    
    echo "Staging $SRR_ID ($PROJECT_ID) ..."
    
    sqlite3 "$DB_FILE" "UPDATE srr_queue SET status='submitted', run_status='SLURM QUEUED', last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR_ID';"
    
    # Slurm pipeline explicitly passing the local target efficiently securely physically
    JOB_ID=$(sbatch --parsable "$DIR/v2_master_node.sh" "$PROJECT_ID" "$SRR_ID" "PAIRED" "$FASTQ_PATH")
    
    if [ -z "$JOB_ID" ]; then
        echo "Failed to submit $SRR_ID. Executing native SQLite rollback..."
        sqlite3 "$DB_FILE" "UPDATE srr_queue SET status='todo', run_status='FETCHED', failure_count=failure_count+1 WHERE srr_id='$SRR_ID';"
    else
        echo "Successfully allocated array $SRR_ID -> Handed securely to SLURM ID: $JOB_ID"
        sqlite3 "$DB_FILE" "UPDATE srr_queue SET slurm_job_id='$JOB_ID' WHERE srr_id='$SRR_ID';"
    fi
done <<< "$JOBS"
