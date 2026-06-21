#!/bin/bash
# fetch_ena_aspera.sh — Download FASTQ from ENA via Aspera (high-speed)
# Usage: ./fetch_ena_aspera.sh <SRR_ID> <TARGET_DIR>

SRR_ID=$1
TARGET_DIR=$2

if [ -z "$SRR_ID" ] || [ -z "$TARGET_DIR" ]; then
    echo "Usage: ./fetch_ena_aspera.sh <SRR_ID> <TARGET_DIR>" >&2
    exit 1
fi

ASCP_KEY="/scratch/hpc-prf-omiks/ja/miniconda3/etc/asperaweb_id_dsa.openssh"

# Get FTP paths from ENA API and convert to Aspera paths
ENA_API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR_ID}&result=read_run&fields=fastq_ftp"
FTP_PATHS=$(curl -s "$ENA_API_URL" | awk 'NR==2 {print $2}')

if [ -z "$FTP_PATHS" ]; then
    echo "ERROR: No ENA paths found for $SRR_ID" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR"

FASTQ_FILES=()
FAIL=0

IFS=';' read -ra URL_ARR <<< "$FTP_PATHS"
for url in "${URL_ARR[@]}"; do
    if [ -n "$url" ]; then
        filename="${url##*/}"
        if [ -s "$TARGET_DIR/$filename" ] && gzip -t "$TARGET_DIR/$filename" 2>/dev/null; then
            echo "Reusing existing $filename" >&2
            FASTQ_FILES+=("$TARGET_DIR/$filename")
            continue
        fi

        # Convert FTP path to Aspera path:
        # ftp.sra.ebi.ac.uk/vol1/fastq/... -> era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/...
        ASPERA_PATH=$(echo "$url" | sed 's|ftp.sra.ebi.ac.uk/|era-fasp@fasp.sra.ebi.ac.uk:/|')

        echo "Downloading via Aspera: $ASPERA_PATH" >&2
        ascp -QT -k1 -l 300m -P 33001 \
            -i "$ASCP_KEY" \
            "$ASPERA_PATH" \
            "$TARGET_DIR/" >&2

        if [ $? -ne 0 ]; then
            echo "ERROR: Aspera download failed for $filename" >&2
            FAIL=1
            break
        fi

        FASTQ_FILES+=("$TARGET_DIR/$filename")
    fi
done

if [ "$FAIL" -eq 1 ] || [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    rm -f "${FASTQ_FILES[@]}" 2>/dev/null
    exit 1
fi

echo "SUCCESS: Aspera download completed for $SRR_ID" >&2
echo "${FASTQ_FILES[@]}"
exit 0
