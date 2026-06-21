#!/bin/bash
# v2_trim.sh
# Usage: ./v2_trim.sh <MODE> <SRR_ID> <LAYOUT> <PHRED> <TARGET_DIR> <FASTQ_PATH1> [FASTQ_PATH2]

set -o pipefail

MODE=$1
SRR_ID=$2
LAYOUT=$3
PHRED=$4
TARGET_DIR=$5
FASTQ1=$6
FASTQ2=$7

BASE="/scratch/hpc-prf-omiks/ja/omiks_project"
ADAPTER="$BASE/resources/adapter_sequences.fasta"


mkdir -p "$TARGET_DIR"
export _JAVA_OPTIONS="-Xmx100g"
TRIM_LOG="$TARGET_DIR/${SRR_ID}_${MODE}_trimmomatic.log"

append_trim_stats() {
    local log_file=$1
    local stats_file=${TRIM_STATS_FILE:-}

    [ -n "$stats_file" ] || return 0
    [ -s "$log_file" ] || return 0

    mkdir -p "$(dirname "$stats_file")"

    python3 - "$log_file" "$stats_file" "${PROJECT_ID:-unknown}" "$SRR_ID" "$MODE" "$LAYOUT" <<'PY'
import re
import sys
from pathlib import Path

log_path, stats_path, project_id, srr_id, mode, layout = sys.argv[1:]
text = Path(log_path).read_text(errors="replace")
stats_file = Path(stats_path)

    header = "project_id\tsrr_id\tmode\tlayout\tinput_reads\tsurviving\tsurviving_pct\tdropped\tdropped_pct\n"
    existing = []
    if stats_file.exists():
        existing = stats_file.read_text(encoding="utf-8", errors="replace").splitlines()
    else:
        existing = [header.rstrip("\n")]

row = None
if layout == "PE":
    m = re.search(
        r"Input Read Pairs:\s+(\d+)\s+Both Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+"
        r"Forward Only Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+"
        r"Reverse Only Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Dropped:\s+(\d+)\s+\(([\d.]+)%\)",
        text,
    )
    if m:
        input_reads = int(m.group(1))
        both = int(m.group(2))
        forward = int(m.group(4))
        reverse = int(m.group(6))
        surviving = both + forward + reverse
        dropped = int(m.group(8))
        surviving_pct = f"{(surviving / input_reads * 100) if input_reads else 0:.2f}"
        dropped_pct = m.group(9)
        row = [
            project_id,
            srr_id,
            mode,
            layout,
            str(input_reads),
            str(surviving),
            surviving_pct,
            str(dropped),
            dropped_pct,
        ]
else:
    m = re.search(
        r"Input Reads:\s+(\d+)\s+Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Dropped:\s+(\d+)\s+\(([\d.]+)%\)",
        text,
    )
    if m:
        row = [project_id, srr_id, mode, layout, *m.groups()]

if row:
    filtered = [existing[0]]
    for line in existing[1:]:
        parts = line.split("\t")
        if len(parts) >= 3 and not (parts[0] == project_id and parts[1] == srr_id and parts[2] == mode):
            filtered.append(line)
    filtered.append("\t".join(row))
    stats_file.write_text("\n".join(filtered) + "\n", encoding="utf-8")
PY
}

if [ "$MODE" = "adapter_only" ]; then
    TRIM_ARGS="ILLUMINACLIP:$ADAPTER:2:30:10"
else
    # Automatically extracts formatting (P5 -> 5)
    P_VAL=${MODE#P}
    TRIM_ARGS="ILLUMINACLIP:$ADAPTER:2:30:10 SLIDINGWINDOW:4:$P_VAL LEADING:$P_VAL TRAILING:$P_VAL"
fi

FASTQ_FILES=()

# Limit Java heap heavily: Because these modes run sequentially, restricting memory
# strictly up to 64GB operates perfectly within the single node 192GB memory constraint!
export JAVA_TOOL_OPTIONS="-Xmx140g"
echo "Executing Trimmomatic in $MODE mode..." >&2

if [ "$LAYOUT" = "PE" ]; then
    OUT_1_PAIRED="$TARGET_DIR/${SRR_ID}_${MODE}_1_paired.fastq.gz"
    OUT_1_UNPAIRED="$TARGET_DIR/${SRR_ID}_${MODE}_1_unpaired.fastq.gz"
    OUT_2_PAIRED="$TARGET_DIR/${SRR_ID}_${MODE}_2_paired.fastq.gz"
    OUT_2_UNPAIRED="$TARGET_DIR/${SRR_ID}_${MODE}_2_unpaired.fastq.gz"
    
    trimmomatic PE -threads ${ACTIVE_THREADS:-12} -$PHRED \
        "$FASTQ1" "$FASTQ2" \
        "$OUT_1_PAIRED" "$OUT_1_UNPAIRED" \
        "$OUT_2_PAIRED" "$OUT_2_UNPAIRED" \
        $TRIM_ARGS 2>&1 | tee "$TRIM_LOG" >&2
        
    if [ $? -ne 0 ]; then echo "ERROR: Trimmomatic string failed execution" >&2; exit 1; fi
    FASTQ_FILES=("$OUT_1_PAIRED" "$OUT_2_PAIRED")
    
    # Rapid cleanup to suppress RAM overflow targets from unpaired sequences
    rm -f "$OUT_1_UNPAIRED" "$OUT_2_UNPAIRED"
else
    OUT_SE="$TARGET_DIR/${SRR_ID}_${MODE}.fastq.gz"
    trimmomatic SE -threads ${ACTIVE_THREADS:-12} -$PHRED \
        "$FASTQ1" \
        "$OUT_SE" \
        $TRIM_ARGS 2>&1 | tee "$TRIM_LOG" >&2
        
    if [ $? -ne 0 ]; then echo "ERROR: Trimmomatic string failed execution" >&2; exit 1; fi
    FASTQ_FILES=("$OUT_SE")
fi

append_trim_stats "$TRIM_LOG"

# Secure path pass back via array intercept standard output
echo "${FASTQ_FILES[@]}"
exit 0
