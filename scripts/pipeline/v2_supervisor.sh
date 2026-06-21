#!/bin/bash
# v2_supervisor.sh — Maintains up to 3 concurrent fetch jobs per endpoint.
# Pipeline jobs are unlimited.
# Usage: screen -dmS supervisor bash v2_supervisor.sh

DB_FILE="$HOME/srr_queue.db"
DIR="/pc2/users/o/omiks001/scripts"
mkdir -p "$DIR/fetch_logs"
MAX_FETCH_FAILURES=10

ENDPOINTS=("ena_aspera")
MAX_PER_EP=3

claim_job() {
    local ep="$1"
    local row
    row=$(sqlite3 -cmd ".timeout 30000" "$DB_FILE" "
        SELECT s.project_id || '|' || s.srr_id
        FROM srr_queue s
        WHERE s.status='todo'
          AND COALESCE(s.failure_count,0) < $MAX_FETCH_FAILURES
          AND (s.run_status IS NULL OR s.run_status='todo' OR s.run_status='FETCH_FAILED' OR s.run_status='INITIALIZED')
          AND NOT EXISTS (
            SELECT 1 FROM fetch_attempts fa_running
            WHERE fa_running.srr_id = s.srr_id
              AND fa_running.status='running'
          )
          AND s.srr_id NOT IN (
            SELECT fa.srr_id FROM fetch_attempts fa
            WHERE fa.endpoint='$ep' AND fa.status='failed'
          )
        ORDER BY s.failure_count ASC, s.last_updated ASC, s.srr_id ASC
        LIMIT 1;")

    [ -z "$row" ] && return 1

    local prj srr changed
    prj=$(echo "$row" | cut -d '|' -f 1)
    srr=$(echo "$row" | cut -d '|' -f 2)

    changed=$(sqlite3 -cmd ".timeout 30000" "$DB_FILE" "
        UPDATE srr_queue
        SET run_status='FETCHING_FILES',
            fetch_lane='$ep',
            last_updated=CURRENT_TIMESTAMP
        WHERE srr_id='$srr'
          AND status='todo'
          AND COALESCE(failure_count,0) < $MAX_FETCH_FAILURES
          AND (run_status IS NULL OR run_status='todo' OR run_status='FETCH_FAILED' OR run_status='INITIALIZED');
        SELECT changes();")

    [ "$changed" -eq 1 ] || return 1
    echo "$prj|$srr"
}

echo "=========================================================="
echo "V2 Supervisor started at $(date)"
echo "Endpoints: ${ENDPOINTS[*]} (max $MAX_PER_EP each)"
echo "=========================================================="

while true; do
    python3 "$DIR/v2_reconcile_queue.py" --db "$DB_FILE" --apply >/dev/null 2>&1 || true

    TODO=$(sqlite3 -cmd ".timeout 30000" "$DB_FILE" \
        "SELECT COUNT(*) FROM srr_queue
         WHERE status='todo'
           AND COALESCE(failure_count,0) < $MAX_FETCH_FAILURES
           AND (run_status IS NULL OR run_status='todo' OR run_status='FETCH_FAILED' OR run_status='INITIALIZED');")

    TOTAL_FETCH=$(squeue -u "$USER" -n v2f_ena_aspera,v2f_ena_http,v2f_ncbi,v2f_aws,v2f_gcp -h -t R,PD 2>/dev/null | wc -l)

    echo "[$(date)] Fetch jobs in SLURM: $TOTAL_FETCH | Queue remaining: $TODO"

    if [ "$TODO" -eq 0 ]; then
        echo "Queue empty. Sleeping 5 minutes."
        sleep 300
        continue
    fi

    SUBMITTED=0

    for EP in "${ENDPOINTS[@]}"; do
        # Count SLURM jobs for this specific endpoint by grepping fetch logs
        # We use a per-endpoint job name so squeue can filter
        ACTIVE=$(squeue -u "$USER" -n "v2f_${EP}" -h -t R,PD 2>/dev/null | wc -l)
        SLOTS=$((MAX_PER_EP - ACTIVE))

        if [ "$SLOTS" -le 0 ]; then
            echo "  [$EP] full ($ACTIVE/$MAX_PER_EP)"
            continue
        fi

        echo "  [$EP] $ACTIVE/$MAX_PER_EP active — filling $SLOTS slots"

        for ((i=0; i<SLOTS; i++)); do
            JOB=$(claim_job "$EP")

            if [ -z "$JOB" ]; then
                echo "    [$EP] no more jobs in queue"
                break
            fi

            PRJ=$(echo "$JOB" | cut -d '|' -f 1)
            SRR=$(echo "$JOB" | cut -d '|' -f 2)

            echo "    [$EP] submitting $SRR ($PRJ)"
            FETCH_JOB_ID=$(sbatch --parsable --job-name="v2f_${EP}" "$DIR/v2_slurm_fetch.sh" "$PRJ" "$SRR" "$EP")

            if [ -z "$FETCH_JOB_ID" ]; then
                echo "    [$EP] submit failed for $SRR; resetting row"
                sqlite3 -cmd ".timeout 30000" "$DB_FILE" "
                    UPDATE srr_queue
                    SET run_status='FETCH_FAILED',
                        fetch_lane=NULL,
                        slurm_job_id=NULL,
                        last_updated=CURRENT_TIMESTAMP
                    WHERE srr_id='$SRR';"
                continue
            fi

            sqlite3 -cmd ".timeout 30000" "$DB_FILE" "
                UPDATE srr_queue
                SET slurm_job_id='$FETCH_JOB_ID',
                    last_updated=CURRENT_TIMESTAMP
                WHERE srr_id='$SRR';"
            SUBMITTED=$((SUBMITTED + 1))
        done
    done

    echo "  Submitted $SUBMITTED new jobs this cycle"
    sleep 120
done
