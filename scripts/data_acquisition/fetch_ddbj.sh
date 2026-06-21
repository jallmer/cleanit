#!/bin/bash
# fetch_ddbj.sh
# Usage: ./fetch_ddbj.sh <SRR_ID> <TARGET_DIR>

SRR_ID=$1
TARGET_DIR=$2

if [ -z "$SRR_ID" ] || [ -z "$TARGET_DIR" ]; then
    echo "Usage: ./fetch_ddbj.sh <SRR_ID> <TARGET_DIR>" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR"

SRR_PREFIX=${SRR_ID:0:6}
ID_LEN=${#SRR_ID}

if [ "$ID_LEN" -eq 9 ]; then
    BASE_URL="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/${SRR_PREFIX}/${SRR_ID}"
elif [ "$ID_LEN" -eq 10 ]; then
    LAST_DIGIT=${SRR_ID:9:1}
    BASE_URL="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/${SRR_PREFIX}/00${LAST_DIGIT}/${SRR_ID}"
elif [ "$ID_LEN" -eq 11 ]; then
    LAST_DIGITS=${SRR_ID:9:2}
    BASE_URL="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/${SRR_PREFIX}/0${LAST_DIGITS}/${SRR_ID}"
else
    echo "ERROR: SRR ID length structure uniquely unsupported by DDBJ mapping block." >&2
    exit 1
fi

FASTQ_FILES=()

# Probing headers mathematically bypasses missing file blocks
wget -q --spider "$BASE_URL/${SRR_ID}_1.fastq.bz2" >&2
if [ $? -eq 0 ]; then
    LAYOUT="PE"
else
    wget -q --spider "$BASE_URL/${SRR_ID}.fastq.bz2" >&2
    if [ $? -eq 0 ]; then
        LAYOUT="SE"
    else
        echo "ERROR: DDBJ DRA repository unavailable." >&2
        exit 1
    fi
fi

if [ "$LAYOUT" = "PE" ]; then
    echo "Acquiring DDBJ BZ2 Array Paths in PARALLEL..." >&2
    wget -q --show-progress "$BASE_URL/${SRR_ID}_1.fastq.bz2" -O "$TARGET_DIR/${SRR_ID}_1.fastq.bz2" >&2 &
    PID1=$!
    
    wget -q --show-progress "$BASE_URL/${SRR_ID}_2.fastq.bz2" -O "$TARGET_DIR/${SRR_ID}_2.fastq.bz2" >&2 &
    PID2=$!
    
    wait $PID1 || exit 1
    wait $PID2 || exit 1
    
    echo "Transcoding DDBJ bz2 vectors securely bounding CPU multithreading natively inside RAM layer..." >&2
    (bzcat "$TARGET_DIR/${SRR_ID}_1.fastq.bz2" | gzip -c > "$TARGET_DIR/${SRR_ID}_1.fastq.gz" && rm -f "$TARGET_DIR/${SRR_ID}_1.fastq.bz2") &
    TPID1=$!
    
    (bzcat "$TARGET_DIR/${SRR_ID}_2.fastq.bz2" | gzip -c > "$TARGET_DIR/${SRR_ID}_2.fastq.gz" && rm -f "$TARGET_DIR/${SRR_ID}_2.fastq.bz2") &
    TPID2=$!
    
    wait $TPID1 || exit 1
    wait $TPID2 || exit 1
    
    FASTQ_FILES=("$TARGET_DIR/${SRR_ID}_1.fastq.gz" "$TARGET_DIR/${SRR_ID}_2.fastq.gz")
else
    echo "Acquiring DDBJ BZ2 Array SE..." >&2
    wget -q --show-progress "$BASE_URL/${SRR_ID}.fastq.bz2" -O "$TARGET_DIR/${SRR_ID}.fastq.bz2" >&2 || exit 1
    
    echo "Transcoding DDBJ bz2 sequence natively into compatible gz formats..." >&2
    bzcat "$TARGET_DIR/${SRR_ID}.fastq.bz2" | gzip -c > "$TARGET_DIR/${SRR_ID}.fastq.gz"
    
    rm -f "$TARGET_DIR/${SRR_ID}.fastq.bz2"
    FASTQ_FILES=("$TARGET_DIR/${SRR_ID}.fastq.gz")
fi

echo "SUCCESS: DDBJ explicitly transposed arrays parallel mapped successfully for $SRR_ID" >&2
echo "${FASTQ_FILES[@]}"
exit 0
