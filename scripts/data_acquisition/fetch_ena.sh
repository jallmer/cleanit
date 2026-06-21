#!/bin/bash
# fetch_ena.sh
# Usage: ./fetch_ena.sh <SRR_ID> <TARGET_DIR>

SRR_ID=$1
TARGET_DIR=$2

if [ -z "$SRR_ID" ] || [ -z "$TARGET_DIR" ]; then
    echo "Usage: ./fetch_ena.sh <SRR_ID> <TARGET_DIR>" >&2
    exit 1
fi

ENA_API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR_ID}&result=read_run&fields=fastq_ftp"
FTP_PATHS=$(curl -s "$ENA_API_URL" | awk 'NR==2 {print $2}')

if [ -z "$FTP_PATHS" ]; then
    echo "ERROR: Failed to retrieve FTP paths for $SRR_ID from ENA API." >&2
    exit 1
fi

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR" || exit 1

FASTQ_FILES=()
PIDS=()

IFS=';' read -ra URL_ARR <<< "$FTP_PATHS"
for url in "${URL_ARR[@]}"; do
    if [ -n "$url" ]; then
        filename="${url##*/}"
        if [ -s "$filename" ] && gzip -t "$filename" 2>/dev/null; then
            echo "Reusing existing $filename" >&2
            FASTQ_FILES+=("$TARGET_DIR/$filename")
            continue
        fi

        rm -f "$filename"
        echo "Downloading ftp://$url to $filename in native PARALLEL..." >&2
        
        timeout 18h wget -q --show-progress --timeout=60 --tries=3 "ftp://$url" -O "$filename" >&2 &
        PIDS+=($!)
        
        FASTQ_FILES+=("$TARGET_DIR/$filename")
    fi
done

# Mathematically block node progression until all parallel networking layers successfully evaluate
for pid in "${PIDS[@]}"; do
    wait $pid
    if [ $? -ne 0 ]; then
        echo "ERROR: Parallel wget stream failed" >&2
        exit 1
    fi
done

echo "SUCCESS: ENA parallel arrays structurally closed for $SRR_ID" >&2
echo "${FASTQ_FILES[@]}"
exit 0
