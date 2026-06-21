#!/bin/bash
# fetch_ncbi.sh — Download only. No extraction. Pipeline handles conversion.
# Usage: ./fetch_ncbi.sh <SRR_ID> <TARGET_DIR>

SRR_ID=$1
TARGET_DIR=$2

if [ -z "$SRR_ID" ] || [ -z "$TARGET_DIR" ]; then
    echo "Usage: ./fetch_ncbi.sh <SRR_ID> <TARGET_DIR>" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR" || exit 1

if [ -s "$TARGET_DIR/${SRR_ID}.sra" ]; then
    echo "Reusing existing ${SRR_ID}.sra" >&2
    echo "$TARGET_DIR/${SRR_ID}.sra"
    exit 0
fi

echo "Prefetching NCBI SRA for $SRR_ID ..." >&2
timeout 12h prefetch --max-size 100G --output-directory "$TARGET_DIR/" "$SRR_ID" >&2

SRA_FILE="$TARGET_DIR/${SRR_ID}/${SRR_ID}.sra"
if [ $? -ne 0 ] || [ ! -f "$SRA_FILE" ]; then
    echo "ERROR: prefetch failed for $SRR_ID" >&2
    exit 1
fi

# Move .sra to target dir root for cleaner paths
mv "$SRA_FILE" "$TARGET_DIR/${SRR_ID}.sra"
rmdir "$TARGET_DIR/${SRR_ID}" 2>/dev/null

echo "SUCCESS: Downloaded $SRR_ID.sra ($(du -sh "$TARGET_DIR/${SRR_ID}.sra" | cut -f1))" >&2

# Output the .sra path (pipeline will handle extraction)
echo "$TARGET_DIR/${SRR_ID}.sra"
exit 0
