#!/bin/bash
# fetch_ncbi_cloud.sh — Download .sra from NCBI cloud (S3/GS). No extraction.
# Usage: ./fetch_ncbi_cloud.sh <SRR_ID> <TARGET_DIR> <CLOUD_SERVICE (s3|gs)>

SRR_ID=$1
TARGET_DIR=$2
CLOUD_SERVICE=${3:-s3}

if [ -z "$SRR_ID" ] || [ -z "$TARGET_DIR" ]; then
    echo "Usage: ./fetch_ncbi_cloud.sh <SRR_ID> <TARGET_DIR> [s3|gs]" >&2
    exit 1
fi

echo "Probing NCBI Locator API for $CLOUD_SERVICE link on $SRR_ID..." >&2

URL=$(curl -s "https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc=${SRR_ID}" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    for res in data.get('result', []):
        for f in res.get('files', []):
            for loc in f.get('locations', []):
                if loc.get('service') == '$CLOUD_SERVICE':
                    print(loc.get('link'))
                    sys.exit(0)
except Exception:
    pass
")

if [ -z "$URL" ]; then
    echo "ERROR: No $CLOUD_SERVICE URL found for $SRR_ID" >&2
    exit 1
fi

echo "Downloading $SRR_ID from $CLOUD_SERVICE..." >&2
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR" || exit 1

if [ -s "${SRR_ID}.sra" ]; then
    echo "Reusing existing ${SRR_ID}.sra" >&2
    echo "$TARGET_DIR/${SRR_ID}.sra"
    exit 0
fi

timeout 12h wget -q --show-progress "$URL" -O "${SRR_ID}.sra" >&2

if [ $? -ne 0 ] || [ ! -f "${SRR_ID}.sra" ] || [ ! -s "${SRR_ID}.sra" ]; then
    echo "ERROR: Download failed for $SRR_ID via $CLOUD_SERVICE" >&2
    rm -f "${SRR_ID}.sra"
    exit 1
fi

echo "SUCCESS: Downloaded ${SRR_ID}.sra ($(du -sh "${SRR_ID}.sra" | cut -f1)) from $CLOUD_SERVICE" >&2

# Output the .sra path (pipeline will handle extraction)
echo "$TARGET_DIR/${SRR_ID}.sra"
exit 0
