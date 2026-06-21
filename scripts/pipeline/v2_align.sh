#!/bin/bash
# v2_align.sh
# Usage: ./v2_align.sh <SRR_ID> <LAYOUT> <PHRED> <INDEX> <TARGET_DIR> <FASTQ_PATH1> [FASTQ_PATH2]

set -o pipefail

SRR_ID=$1
LAYOUT=$2
PHRED=$3
INDEX=$4
TARGET_DIR=$5
FASTQ1=$6
FASTQ2=$7


mkdir -p "$TARGET_DIR"
BAM_OUT="$TARGET_DIR/${SRR_ID}_sorted.bam"
BOWTIE_LOG="$TARGET_DIR/${SRR_ID}_bowtie2.log"

append_bowtie_stats() {
    local log_file=$1
    local stats_file=${BOWTIE_STATS_FILE:-}

    [ -n "$stats_file" ] || return 0
    [ -s "$log_file" ] || return 0

    mkdir -p "$(dirname "$stats_file")"

    python3 - "$log_file" "$stats_file" "${PROJECT_ID:-unknown}" "$SRR_ID" "${MODE_NAME:-unknown}" "$LAYOUT" <<'PY'
import re
import sys
from pathlib import Path

log_path, stats_path, project_id, srr_id, mode_name, layout = sys.argv[1:]
text = Path(log_path).read_text(errors="replace")
stats_file = Path(stats_path)

header = (
    "project_id\tsrr_id\tmode\tlayout\ttotal_reads\tpaired_reads\tunpaired_reads\taligned_exactly_1\taligned_gt1\taligned_0\t"
    "concordant_exactly_1\tconcordant_gt1\tdiscordant_1\tpairs_0_concordant\toverall_alignment_rate\n"
)
if stats_file.exists():
    existing = stats_file.read_text(encoding="utf-8", errors="replace").splitlines()
else:
    existing = [header.rstrip("\n")]

def grab(pattern):
    m = re.search(pattern, text, re.MULTILINE)
    return m.group(1) if m else ""

row = [project_id, srr_id, mode_name, layout, grab(r"(\d+)\s+reads; of these:")]

if layout == "PE":
    row.extend(
        [
            grab(r"(\d+)\s+\([^)]+\)\s+were paired; of these:"),
            grab(r"(\d+)\s+\([^)]+\)\s+were unpaired; of these:"),
            "",
            "",
            "",
            grab(r"(\d+)\s+\([^)]+\)\s+aligned concordantly exactly 1 time"),
            grab(r"(\d+)\s+\([^)]+\)\s+aligned concordantly >1 times"),
            grab(r"(\d+)\s+\([^)]+\)\s+aligned discordantly 1 time"),
            grab(r"(\d+)(?:\s+\([^)]+\))?\s+pairs aligned concordantly 0 times"),
            grab(r"([\d.]+%)\s+overall alignment rate"),
        ]
    )
else:
    row.extend(
        [
            "",
            "",
            grab(r"(\d+)\s+\([^)]+\)\s+aligned exactly 1 time$"),
            grab(r"(\d+)\s+\([^)]+\)\s+aligned >1 times$"),
            grab(r"(\d+)\s+\([^)]+\)\s+aligned 0 times$"),
            "",
            "",
            "",
            "",
            grab(r"([\d.]+%)\s+overall alignment rate"),
        ]
    )

if row[4] and row[-1]:
    filtered = [existing[0]]
    for line in existing[1:]:
        parts = line.split("\t")
        if len(parts) >= 3 and not (parts[0] == project_id and parts[1] == srr_id and parts[2] == mode_name):
            filtered.append(line)
    filtered.append("\t".join(row))
    stats_file.write_text("\n".join(filtered) + "\n", encoding="utf-8")
PY
}

echo "Executing Bowtie2 Matrix Arrays Mapping to Samtools Index..." >&2

if [ "$LAYOUT" = "PE" ]; then
    bowtie2 -p ${ACTIVE_THREADS:-12} --$PHRED -x "$INDEX" -1 "$FASTQ1" -2 "$FASTQ2" 2> >(tee "$BOWTIE_LOG" >&2) | \
    /scratch/hpc-prf-omiks/fb/omiks_project/envs/samtools_env/bin/samtools view -bS - 2>&2 | \
    /scratch/hpc-prf-omiks/fb/omiks_project/envs/samtools_env/bin/samtools sort -@ ${ACTIVE_THREADS:-12} -T "$TARGET_DIR/tmp_sort_$$" -o "$BAM_OUT" 2>&2
else
    bowtie2 -p ${ACTIVE_THREADS:-12} --$PHRED -x "$INDEX" -U "$FASTQ1" 2> >(tee "$BOWTIE_LOG" >&2) | \
    /scratch/hpc-prf-omiks/fb/omiks_project/envs/samtools_env/bin/samtools view -bS - 2>&2 | \
    /scratch/hpc-prf-omiks/fb/omiks_project/envs/samtools_env/bin/samtools sort -@ ${ACTIVE_THREADS:-12} -T "$TARGET_DIR/tmp_sort_$$" -o "$BAM_OUT" 2>&2
fi

if [ $? -ne 0 ] || [ ! -f "$BAM_OUT" ]; then
    echo "ERROR: Bowtie2 to Samtools mapping stream failed execution" >&2
    exit 1
fi

append_bowtie_stats "$BOWTIE_LOG"

# Structurally pass variable over array intercept standard output
echo "$BAM_OUT"
exit 0
