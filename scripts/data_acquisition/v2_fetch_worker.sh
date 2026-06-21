#!/bin/bash
# v2_fetch_worker.sh
# Usage: bash v2_fetch_worker.sh <PROJECT_ID> <SRR_ID> <ENDPOINT>

PROJECT=$1
SRR=$2
ENDPOINT=$3
DB_FILE="$HOME/srr_queue.db"
DIR="/pc2/users/o/omiks001/scripts"
MAX_FETCH_FAILURES=10

sq() { sqlite3 -cmd ".timeout 30000" "$@"; }

echo "=========================================================="
echo "Worker [$ENDPOINT] processing: $SRR ($PROJECT)"
echo "=========================================================="

claim_fetch_attempt() {
    sq "$DB_FILE" "
        BEGIN IMMEDIATE;
        INSERT INTO fetch_attempts (srr_id, endpoint, status)
        SELECT '$SRR', '$ENDPOINT', 'running'
        WHERE NOT EXISTS (
            SELECT 1 FROM fetch_attempts
            WHERE srr_id='$SRR' AND status='running'
        );
        SELECT changes();
        COMMIT;" | tail -n 1
}

get_attempt_id() {
    sq "$DB_FILE" "
        SELECT id
        FROM fetch_attempts
        WHERE srr_id='$SRR'
          AND endpoint='$ENDPOINT'
          AND status='running'
        ORDER BY id DESC
        LIMIT 1;"
}

# Acquire the fetch slot atomically. If another worker already owns this SRR,
# exit before touching the shared target directory.
CLAIMED=$(claim_fetch_attempt)

if [ "$CLAIMED" != "1" ]; then
    echo "[$ENDPOINT] Duplicate fetch prevented for $SRR; another fetch attempt is already running"
    exit 0
fi

ATTEMPT_ID=$(get_attempt_id)

if [ -z "$ATTEMPT_ID" ]; then
    echo "[$ENDPOINT] ERROR: fetch attempt claim succeeded but attempt id lookup failed for $SRR"
    sq "$DB_FILE" "UPDATE srr_queue SET run_status='FETCH_FAILED', fetch_lane=NULL, slurm_job_id=NULL, last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"
    exit 1
fi

TARGET_DIR="/scratch/hpc-prf-omiks/ja/temp/$PROJECT/$SRR"
mkdir -p "$TARGET_DIR"

sq "$DB_FILE" "UPDATE srr_queue SET run_status='FETCHING_FILES', fetch_lane='$ENDPOINT', last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"

# --- FETCH ---
FASTQ_P=$(bash "$DIR/fetch_data.sh" "$SRR" "$TARGET_DIR" "$ENDPOINT")

if [ $? -ne 0 ]; then
    echo "[$ENDPOINT] FAILED: download returned non-zero for $SRR"
    sq "$DB_FILE" "UPDATE fetch_attempts SET status='failed', error_msg='download_error', finished_at=CURRENT_TIMESTAMP WHERE id=$ATTEMPT_ID;"
    sq "$DB_FILE" "UPDATE srr_queue
        SET run_status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'FETCH_ABANDONED' ELSE 'FETCH_FAILED' END,
            fetch_lane='$ENDPOINT',
            failure_count=COALESCE(failure_count,0)+1,
            status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'abandoned' ELSE 'todo' END,
            slurm_job_id=NULL,
            last_updated=CURRENT_TIMESTAMP
        WHERE srr_id='$SRR';"
    exit 1
fi

echo "[$ENDPOINT] Download finished for $SRR"

# --- Determine what we got: .sra or .fastq.gz ---
SRA_FILE=$(ls -1 "$TARGET_DIR"/*.sra 2>/dev/null | head -1)
if [ -n "$SRA_FILE" ]; then
    # NCBI path: got .sra file, validate it
    REAL_PATHS="$SRA_FILE"
    echo "[$ENDPOINT] Got .sra file: $SRA_FILE"
    TOTAL_BYTES=$(stat -c%s "$SRA_FILE")
    IS_SRA=1

    if [ "$TOTAL_BYTES" -eq 0 ]; then
        echo "[$ENDPOINT] FATAL: .sra file is 0 bytes"
        rm -f "$SRA_FILE"
        sq "$DB_FILE" "UPDATE fetch_attempts SET status='failed', error_msg='empty_sra', bytes_downloaded=0, finished_at=CURRENT_TIMESTAMP WHERE id=$ATTEMPT_ID;"
        sq "$DB_FILE" "UPDATE srr_queue
            SET run_status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'FETCH_ABANDONED' ELSE 'FETCH_FAILED' END,
                fetch_lane='$ENDPOINT',
                failure_count=COALESCE(failure_count,0)+1,
                status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'abandoned' ELSE 'todo' END,
                slurm_job_id=NULL,
                last_updated=CURRENT_TIMESTAMP
            WHERE srr_id='$SRR';"
        exit 1
    fi

    # Validate .sra integrity
    if command -v vdb-validate &>/dev/null; then
        echo "[$ENDPOINT] Validating .sra with vdb-validate..."
        if ! vdb-validate "$SRA_FILE" 2>&1; then
            echo "[$ENDPOINT] FATAL: .sra failed vdb-validate"
            rm -f "$SRA_FILE"
            sq "$DB_FILE" "UPDATE fetch_attempts SET status='failed', error_msg='sra_corrupt', bytes_downloaded=$TOTAL_BYTES, finished_at=CURRENT_TIMESTAMP WHERE id=$ATTEMPT_ID;"
            sq "$DB_FILE" "UPDATE srr_queue
                SET run_status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'FETCH_ABANDONED' ELSE 'FETCH_FAILED' END,
                    fetch_lane='$ENDPOINT',
                    failure_count=COALESCE(failure_count,0)+1,
                    status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'abandoned' ELSE 'todo' END,
                    slurm_job_id=NULL,
                    last_updated=CURRENT_TIMESTAMP
                WHERE srr_id='$SRR';"
            exit 1
        fi
    fi
else
    # ENA/cloud path: got .fastq.gz files
    IS_SRA=0
    # Gzip any raw .fastq first
    for f in "$TARGET_DIR"/*.fastq; do
        [ -f "$f" ] && gzip "$f"
    done

    REAL_PATHS=$(ls -1d "$TARGET_DIR"/*.fastq.gz 2>/dev/null | tr '\n' ' ' | xargs)

    if [ -z "$REAL_PATHS" ]; then
        echo "[$ENDPOINT] FATAL: no .fastq.gz or .sra files found"
        sq "$DB_FILE" "UPDATE fetch_attempts SET status='failed', error_msg='no_files', finished_at=CURRENT_TIMESTAMP WHERE id=$ATTEMPT_ID;"
        sq "$DB_FILE" "UPDATE srr_queue
            SET run_status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'FETCH_ABANDONED' ELSE 'FETCH_FAILED' END,
                fetch_lane='$ENDPOINT',
                failure_count=COALESCE(failure_count,0)+1,
                status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'abandoned' ELSE 'todo' END,
                slurm_job_id=NULL,
                last_updated=CURRENT_TIMESTAMP
            WHERE srr_id='$SRR';"
        exit 1
    fi

    # Validate .fastq.gz files
    echo "[$ENDPOINT] Validating downloaded files..."
    VALID=1
    FAIL_REASON=""

    ENA_BYTES=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_bytes" \
        | awk 'NR==2 {print $2}')

    if [ -n "$ENA_BYTES" ]; then
        IFS=';' read -ra EXPECTED <<< "$ENA_BYTES"
        EXPECTED_SORTED=$(printf '%s\n' "${EXPECTED[@]}" | sort -n)
        ACTUAL_SORTED=$(for f in $REAL_PATHS; do stat -c%s "$f"; done | sort -n)
        if [ "$EXPECTED_SORTED" != "$ACTUAL_SORTED" ]; then
            VALID=0
            FAIL_REASON="size_mismatch"
        fi
    fi

    # Always check gzip integrity regardless of size match
    if [ "$VALID" -eq 1 ]; then
        for f in $REAL_PATHS; do
            if ! gzip -t "$f" 2>/dev/null; then
                VALID=0
                FAIL_REASON="gzip_corrupt"
                break
            fi
        done
    fi

    for f in $REAL_PATHS; do
        if [ "$(stat -c%s "$f")" -eq 0 ]; then
            VALID=0
            FAIL_REASON="empty_file"
        fi
    done

    TOTAL_BYTES=$(du -cb $REAL_PATHS | grep total | awk '{print $1}')

    if [ "$VALID" -eq 0 ]; then
        echo "[$ENDPOINT] Validation FAILED: $FAIL_REASON"
        rm -f $REAL_PATHS
        sq "$DB_FILE" "UPDATE fetch_attempts SET status='failed', error_msg='$FAIL_REASON', bytes_downloaded=$TOTAL_BYTES, finished_at=CURRENT_TIMESTAMP WHERE id=$ATTEMPT_ID;"
        sq "$DB_FILE" "UPDATE srr_queue
            SET run_status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'FETCH_ABANDONED' ELSE 'FETCH_FAILED' END,
                fetch_lane='$ENDPOINT',
                failure_count=COALESCE(failure_count,0)+1,
                status=CASE WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'abandoned' ELSE 'todo' END,
                slurm_job_id=NULL,
                last_updated=CURRENT_TIMESTAMP
            WHERE srr_id='$SRR';"
        exit 1
    fi
    echo "[$ENDPOINT] Validation PASSED"
fi

# --- SUCCESS: log attempt, update queue, submit pipeline ---
sq "$DB_FILE" "UPDATE fetch_attempts SET status='success', bytes_downloaded=$TOTAL_BYTES, finished_at=CURRENT_TIMESTAMP WHERE id=$ATTEMPT_ID;"

SIZE_GB=$(( TOTAL_BYTES / 1024 / 1024 / 1024 + 1 ))
if [ "$IS_SRA" = "1" ]; then
    sq "$DB_FILE" "UPDATE srr_queue SET local_fastq_path='$REAL_PATHS', run_status='FETCHED', fetch_lane='$ENDPOINT', size_gb=$SIZE_GB, sra_size_bytes=$TOTAL_BYTES, last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"
else
    sq "$DB_FILE" "UPDATE srr_queue SET local_fastq_path='$REAL_PATHS', run_status='FETCHED', fetch_lane='$ENDPOINT', size_gb=$SIZE_GB, fastq_size_bytes=$TOTAL_BYTES, last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"
fi
sq "$DB_FILE" "UPDATE srr_queue SET status='submitted', run_status='SLURM QUEUED', last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"

JOB_ID=$(sbatch --parsable "$DIR/v2_master_node.sh" "$PROJECT" "$SRR" "PAIRED" "$REAL_PATHS")

if [ -z "$JOB_ID" ]; then
    echo "[$ENDPOINT] ERROR: sbatch rejected"
    sq "$DB_FILE" "UPDATE srr_queue SET status='todo', run_status='FETCHED', failure_count=failure_count+1 WHERE srr_id='$SRR';"
    exit 1
else
    sq "$DB_FILE" "UPDATE srr_queue SET slurm_job_id='$JOB_ID' WHERE srr_id='$SRR';"
    echo "[$ENDPOINT] DONE: $SRR submitted as SLURM job $JOB_ID"
    exit 0
fi
