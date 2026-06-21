#!/bin/bash
# v2_data_fetcher.sh
# Usage: bash v2_data_fetcher.sh [LIMIT]
# Purely natively handles the 5 Cloud Endpoints (NCBI, EBI Aspera, EBI FTP, AWS, GCP) without draining SLURM Cores!

LIMIT=${1:-5}
DB_FILE="$HOME/srr_queue.db"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "=========================================================="
echo "V2 Global SQLite Data Acquisition Fetcher"
echo "LIMIT: $LIMIT"
echo "=========================================================="

# Select jobs explicitly needing Physical Data Fetching natively
JOBS=$(sqlite3 "$DB_FILE" "SELECT project_id, srr_id FROM srr_queue WHERE status='todo' AND (run_status != 'FETCHED' OR run_status IS NULL) ORDER BY failure_count ASC LIMIT $LIMIT;")

if [ -z "$JOBS" ]; then
    echo "No unfetched sequences detected in the dynamic queue. Complete!"
    exit 0
fi

for job in $JOBS; do
    PRJ=$(echo "$job" | cut -d'|' -f1)
    SRR=$(echo "$job" | cut -d'|' -f2)
    
    echo "----------------------------------------------------------"
    echo "Targeting SRR: $SRR (Project: $PRJ)"
    
    # 1. Ensure Target Temporary storage bound exists explicitly outside SLURM cache!
    TARGET_DIR="/scratch/hpc-prf-omiks/ja/temp/$PRJ/$SRR"
    mkdir -p "$TARGET_DIR"
    
    # 2. Prevent loop duplication cleanly
    sqlite3 "$DB_FILE" "UPDATE srr_queue SET run_status='FETCHING_FILES', last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"
    
    # 3. Fire the identical 5 Endpoint Array script natively mapped
    FASTQ_P=$(bash "$DIR/fetch_data.sh" "$SRR" "$TARGET_DIR" "random")
    
    if [ $? -eq 0 ]; then
        echo "SUCCESS: $SRR cleanly downloaded perfectly!"
        
        # 4. Mathematically map precise absolute strings for database
        # Make sure they are gzipped if ending in .fastq!
        for f in $TARGET_DIR/*.fastq; do
            if [ -f "$f" ]; then
                echo "Gzipping explicit raw extraction $f natively..."
                gzip "$f"
            fi
        done
        
        REAL_PATHS=$(ls -1d $TARGET_DIR/*.fastq.gz 2>/dev/null | tr '\n' ' ' | xargs)
        
        if [ -z "$REAL_PATHS" ]; then
            echo "FATAL: Fetch reported 0 but physically no arrays found in directory!"
            sqlite3 "$DB_FILE" "UPDATE srr_queue SET run_status='FETCH_FAILED', failure_count=failure_count+1 WHERE srr_id='$SRR';"
        else
            echo "PATHS RESOLVED: $REAL_PATHS"
            # Formally calculate total physical array footprint in GB
            SIZE_BYTES=$(du -cb $REAL_PATHS | grep total | awk '{print $1}')
            SIZE_GB=$(( SIZE_BYTES / 1024 / 1024 / 1024 + 1 ))
            
            # Secure execution hand-off perfectly cleanly into the native database!
            sqlite3 "$DB_FILE" "UPDATE srr_queue SET local_fastq_path='$REAL_PATHS', run_status='FETCHED', size_gb=$SIZE_GB, last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"
            echo "DATABASE TRANSACTION EXACTLY SECURED. $SRR is now legally explicitly ready for Submission Arrays!"
        fi
    else
        echo "FAILED: $SRR structurally crashed across all 5 fallback endpoints physically!"
        sqlite3 "$DB_FILE" "UPDATE srr_queue SET run_status='FETCH_FAILED', failure_count=failure_count+1, last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR';"
    fi
    
    echo "----------------------------------------------------------"
done

echo "Fetch Sub-Routine Completed Structurally!"
