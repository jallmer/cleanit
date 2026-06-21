#!/bin/bash
# fetch_ddbj_sra.sh
# Usage: ./fetch_ddbj_sra.sh <SRR_ID> <TARGET_DIR>

SRR_ID=$1
TARGET_DIR=$2

if [ -z "$SRR_ID" ] || [ -z "$TARGET_DIR" ]; then
    echo "Usage: ./fetch_ddbj_sra.sh <SRR_ID> <TARGET_DIR>" >&2
    exit 1
fi

# Obtain the associated SRX Experiment identifier completely autonomously from the ENA database API
ENA_API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR_ID}&result=read_run&fields=experiment_accession"
SRX_ID=$(curl -s "$ENA_API_URL" | awk 'NR==2 {print $2}')

if [ -z "$SRX_ID" ]; then
    echo "ERROR: Failed to computationally map SRX identifier for $SRR_ID natively from ENA metadata API." >&2
    exit 1
fi

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR" || exit 1

SRX_PREFIX=${SRX_ID:0:6}
ID_LEN=${#SRX_ID}

# Complex Japanese Directory Routing string analysis logic bounding explicitly to the SRX parent 
if [ "$ID_LEN" -eq 9 ]; then
    BASE_URL="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/${SRX_PREFIX}/${SRX_ID}/${SRR_ID}"
elif [ "$ID_LEN" -eq 10 ]; then
    LAST_DIGIT=${SRX_ID:9:1}
    BASE_URL="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/${SRX_PREFIX}/00${LAST_DIGIT}/${SRX_ID}/${SRR_ID}"
elif [ "$ID_LEN" -eq 11 ]; then
    LAST_DIGITS=${SRX_ID:9:2}
    BASE_URL="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/${SRX_PREFIX}/0${LAST_DIGITS}/${SRX_ID}/${SRR_ID}"
else
    echo "ERROR: SRX ID mathematical length structure uniquely unsupported by DDBJ mapping block." >&2
    exit 1
fi

echo "Probing DDBJ exclusively for natively mirrored SRA sequential datasets... ($BASE_URL/${SRR_ID}.sralite)" >&2

wget -q --spider "$BASE_URL/${SRR_ID}.sralite" >&2
if [ $? -ne 0 ]; then
    # Systematically check if Japan specifically maintains the raw .sra block outside the litesra parameter
    wget -q --spider "$BASE_URL/${SRR_ID}.sra" >&2
    if [ $? -ne 0 ]; then
        echo "ERROR: DDBJ DRA repository does not structurally host the mapped .sra/.sralite framework for $SRR_ID." >&2
        exit 1
    else
        EXT=".sra"
    fi
else
    EXT=".sralite"
fi

echo "Acquiring massive Japanese DDBJ SRA Mirror sequentially for $SRR_ID ..." >&2
wget -q --show-progress "$BASE_URL/${SRR_ID}${EXT}" -O "${SRR_ID}${EXT}" >&2

if [ $? -ne 0 ] || [ ! -f "${SRR_ID}${EXT}" ]; then
    echo "ERROR: Japanese wget array fundamentally completely failed for $SRR_ID" >&2
    exit 1
fi

echo "Dumping DDBJ Sequence linearly to natively readable fastq variables ($TARGET_DIR) ..." >&2
# If it dynamically evaluates as a .sralite block, forcefully rename it to seamlessly trick fasterq-dump into native processing
if [ "$EXT" = ".sralite" ]; then
    mv "${SRR_ID}.sralite" "${SRR_ID}.sra"
fi


timeout 4h fasterq-dump "${SRR_ID}.sra" -e "${SLURM_CPUS_PER_TASK:-4}" -O "$TARGET_DIR/" -t "$TARGET_DIR/" >&2


FASTQ_FILES=()
num_fastq=$(find "$TARGET_DIR" -type f -name "*.fastq" -size +0c | wc -l)
if [ "$num_fastq" -eq 0 ]; then
    echo "ERROR: Zero physical pipeline yields from DDBJ fasterq-dump extraction matrix for $SRR_ID." >&2
    rm -f "${SRR_ID}.sra"
    exit 1
fi

echo "Mathematically dynamically compressing extracted DDBJ arrays natively into gzip vectors..." >&2
# Map physical parallel computation hooks continuously against the single-node SLURM parameters
PIDS=()
for file in "$TARGET_DIR"/*.fastq; do
    if [ -f "$file" ]; then
        gzip "$file" &
        PIDS+=($!)
        FASTQ_FILES+=("${file}.gz")
    fi
done

for pid in "${PIDS[@]}"; do
    wait $pid
    if [ $? -ne 0 ]; then
        echo "ERROR: DDBJ Parallel GZ Compression block aggressively failed." >&2
        rm -f "${SRR_ID}.sra"
        exit 1
    fi
done

# Violently drop the enormous SRA sequence vector from the Memory disk instantly
rm -f "${SRR_ID}.sra"

echo "SUCCESS: DDBJ securely dumped SRA architecture successfully processing parallel-gzips for $SRR_ID" >&2
echo "${FASTQ_FILES[@]}"
exit 0
