#!/bin/bash
# fetch_data.sh
# Usage: ./fetch_data.sh <SRR_ID> <TARGET_DIR>

SRR_ID=$1
TARGET_DIR=$2
FETCH_MODE=${3:-random}
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "[$(date)] Starting natively targeted data fetch orchestrator for $SRR_ID into $TARGET_DIR (MODE: $FETCH_MODE)" >&2

if [ "$FETCH_MODE" = "ena" ]; then
    SHUFFLED=("fetch_ena.sh")
elif [ "$FETCH_MODE" = "ena_ftp" ]; then
    SHUFFLED=("fetch_ena.sh")
elif [ "$FETCH_MODE" = "ena_http" ]; then
    SHUFFLED=("fetch_ena_http.sh")
elif [ "$FETCH_MODE" = "ena_aspera" ]; then
    SHUFFLED=("fetch_ena_aspera.sh")
elif [ "$FETCH_MODE" = "ncbi" ]; then
    SHUFFLED=("fetch_ncbi.sh")
elif [ "$FETCH_MODE" = "aws" ]; then
    SHUFFLED=("fetch_aws.sh")
elif [ "$FETCH_MODE" = "gcp" ]; then
    SHUFFLED=("fetch_gcp.sh")
else
    echo "Probing global APIs defensively to assemble mathematically valid endpoint array for $SRR_ID..." >&2
    VALID_FETCHERS=()
    
    # 1. Fast Probe ENA API
    ENA_URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR_ID}&result=read_run&fields=fastq_ftp" | awk 'NR==2 {print $2}')
    if [ -n "$ENA_URLS" ]; then
        VALID_FETCHERS+=("fetch_ena_aspera.sh" "fetch_ena_http.sh" "fetch_ena.sh")
    fi
    
    # 2. Deep Probe NCBI JSON API securely via python
    CLOUD_SERVICES=$(curl -s "https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc=${SRR_ID}" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    s = set()
    for res in data.get('result', []):
        for f in res.get('files', []):
            for loc in f.get('locations', []):
                if loc.get('service'): s.add(loc.get('service'))
    print(' '.join(s))
except: pass
")
    
    if [[ "$CLOUD_SERVICES" == *"s3"* ]]; then VALID_FETCHERS+=("fetch_aws.sh"); fi
    if [[ "$CLOUD_SERVICES" == *"gs"* ]]; then VALID_FETCHERS+=("fetch_gcp.sh"); fi
    
    # NCBI SRA is the absolute terminal infrastructure fallback limit natively
    VALID_FETCHERS+=("fetch_ncbi.sh")
    
    echo "Discovered mathematically physically verified endpoints: ${VALID_FETCHERS[*]}" >&2
    
    SHUFFLED=()
    for f in "${VALID_FETCHERS[@]}"; do [ "$f" == "fetch_gcp.sh" ] && SHUFFLED+=("$f"); done
    for f in "${VALID_FETCHERS[@]}"; do [ "$f" == "fetch_aws.sh" ] && SHUFFLED+=("$f"); done
    
    REST=()
    for f in "${VALID_FETCHERS[@]}"; do [ "$f" != "fetch_gcp.sh" ] && [ "$f" != "fetch_aws.sh" ] && REST+=("$f"); done
    
    if [ ${#REST[@]} -gt 0 ]; then
        SHUFFLED+=($(shuf -e "${REST[@]}"))
    fi
fi

echo "[$(date)] Structurally balancing network load: Execution hierarchy randomly defined as [ ${SHUFFLED[*]} ]" >&2

mkdir -p "$TARGET_DIR"

for fetch_script in "${SHUFFLED[@]}"; do
    echo "[$(date)] >>> ATTEMPTING NETWORK BOUNDARY: $fetch_script <<<" >&2

    FASTQ_PATHS=$(bash "$DIR/$fetch_script" "$SRR_ID" "$TARGET_DIR")
    FETCH_EXIT=$?
    
    if [ $FETCH_EXIT -eq 0 ]; then
        echo "[$(date)] SUCCESS: $fetch_script correctly acquired the complete structure." >&2
        echo "$FASTQ_PATHS"
        exit 0
    else
        echo "[$(date)] WARNING: $fetch_script formally failed pipeline execution. Triaging to next immediate fallback parameter..." >&2
    fi
done

echo "[$(date)] FATAL ERROR: Entire shuffled randomization block (DDBJ, ENA FTP, ENA HTTP, NCBI) mathematically failed to yield $SRR_ID." >&2
exit 1
