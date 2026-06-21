#!/bin/bash
# resubmit_pipelines.sh — Resubmit pipeline jobs for entries with files already downloaded
DIR="/pc2/users/o/omiks001/scripts"
DB="/pc2/users/o/omiks001/srr_queue.db"

sqlite3 "$DB" "SELECT srr_id, project_id, local_fastq_path FROM srr_queue WHERE status='submitted' AND run_status='SLURM QUEUED' AND local_fastq_path IS NOT NULL AND local_fastq_path != '';" | while IFS='|' read -r SRR PRJ PATHS; do
  FIRST=$(echo "$PATHS" | awk '{print $1}')
  if [ -f "$FIRST" ]; then
    JOB=$(sbatch --parsable "$DIR/v2_master_node.sh" "$PRJ" "$SRR" "PAIRED" "$PATHS")
    echo "RESUBMITTED: $SRR ($PRJ) -> job $JOB"
  else
    echo "SKIP (files gone): $SRR — resetting to todo"
    sqlite3 "$DB" "UPDATE srr_queue SET status='todo', run_status=NULL, local_fastq_path=NULL WHERE srr_id='$SRR';"
  fi
done
