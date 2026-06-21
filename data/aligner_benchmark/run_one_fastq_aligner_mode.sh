#!/bin/bash
# One task = one raw FASTQ sample x trimming mode x aligner benchmark.
# This runner never calls SRA tools.

#SBATCH -J fastq_align_bench
#SBATCH -o /scratch/hpc-prf-omiks/ja/aligner_benchmark/logs/%x_%A_%a.out
#SBATCH -e /scratch/hpc-prf-omiks/ja/aligner_benchmark/logs/%x_%A_%a.err
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180G
#SBATCH -p normal

set -euo pipefail

CORES="${SLURM_CPUS_PER_TASK:-36}"
TRIM_THREADS="${TRIM_THREADS:-8}"
TASK_ID="${1:-${SLURM_ARRAY_TASK_ID:-1}}"

ROOT="/scratch/hpc-prf-omiks/ja/aligner_benchmark"
TASKS="$ROOT/aligner_fastq_216_tasks.tsv"
RESULTS="$ROOT/results/fastq_aligner_mode_speed.tsv"
PREP_LOG="$ROOT/results/fastq_trim_prep.tsv"
CACHE="$ROOT/fastq_cache"
WORK_ROOT="/scratch/hpc-prf-omiks/ja/tmp/fastq_aligner_benchmark_jobs"

STAR_BIN="/scratch/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/STAR"
STAR_INDEX="/scratch/hpc-prf-omiks/ja/genomes/combined_star_index"
HISAT2_BIN="/scratch/hpc-prf-omiks/ja/envs/hisat2_env/bin/hisat2"
HISAT2_INDEX="/scratch/hpc-prf-omiks/ja/genomes/combined_hisat2_index"
BOWTIE2_BIN="${BOWTIE2_BIN:-/scratch/hpc-prf-omiks/ja/envs/hisat2_env/bin/bowtie2}"
BOWTIE2_INDEX="/scratch/hpc-prf-omiks/mimo/zika_rnaseq_project/reference/human/bowtie2_index/GRCh38_p14"
TRIMMOMATIC="/scratch/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/trimmomatic"
ADAPTER="/scratch/hpc-prf-omiks/ja/omiks_project/resources/adapter_sequences.fasta"

export PATH="/scratch/hpc-prf-omiks/ja/envs/hisat2_env/bin:/scratch/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin:$PATH"
export TMPDIR="$WORK_ROOT/tmp_${SLURM_JOB_ID:-manual}_${TASK_ID}"

mkdir -p "$ROOT/results" "$ROOT/logs" "$CACHE/trimmed_fastq_manifest" "$WORK_ROOT" "$TMPDIR"

init_table() {
    local path=$1
    local header=$2
    if [ ! -s "$path" ]; then
        printf "%s\n" "$header" > "$path"
    fi
}

append_row() {
    local path=$1
    shift
    {
        flock 9
        printf "%s\n" "$*" >> "$path"
    } 9>"$path.lock"
}

init_table "$RESULTS" "timestamp	task_id	project_id	srr_id	condition	cell_line	mode	aligner	layout	total_reads_M	fastq_mb	bowtie_trimmed_sec	bowtie_alignment_rate	threads	duration_sec	exit_code	status	overall_alignment_rate	unique_reads	multimapped_reads	unmapped_reads	log_path"
init_table "$PREP_LOG" "timestamp	task_id	project_id	srr_id	mode	stage	status	threads	duration_sec	exit_code	message"

line=$(awk -v n="$TASK_ID" 'NR == n + 1 {print}' "$TASKS")
if [ -z "$line" ]; then
    echo "No task row for task_id=$TASK_ID" >&2
    exit 1
fi

IFS=$'\t' read -r TASK_ID_FILE PROJECT_ID SRR_ID CONDITION CELL_LINE LAYOUT MODE ALIGNER TOTAL_READS_M BOWTIE_TRIMMED_SEC BOWTIE_ALIGNMENT_RATE FASTQ_MB R1_PATH R2_PATH <<< "$line"
WORK="$WORK_ROOT/${TASK_ID}_${SRR_ID}_${MODE}_${ALIGNER}_${SLURM_JOB_ID:-manual}"
mkdir -p "$WORK"

cleanup() {
    rm -rf "$WORK" "$TMPDIR"
}
trap cleanup EXIT

read_command_for() {
    case "$1" in
        *.gz) printf "zcat" ;;
        *) printf "cat" ;;
    esac
}

find_bowtie2() {
    if [ -n "$BOWTIE2_BIN" ] && [ -x "$BOWTIE2_BIN" ]; then
        printf "%s\n" "$BOWTIE2_BIN"
        return 0
    fi
    if command -v bowtie2 >/dev/null 2>&1; then
        command -v bowtie2
        return 0
    fi
    if [ -s "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        # The original Zika workflow used this environment for Bowtie2.
        # shellcheck disable=SC1091
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
        conda activate zika_rnaseq >/dev/null 2>&1 || true
        if command -v bowtie2 >/dev/null 2>&1; then
            command -v bowtie2
            return 0
        fi
    fi
    return 1
}

ensure_trimmed_fastqs() {
    local mode=$1
    local in1=$2
    local in2=$3
    if [ "$mode" = "untrimmed" ]; then
        printf "%s\t%s\n" "$in1" "$in2"
        return 0
    fi

    local trim_dir="$CACHE/trimmed_fastq_manifest/$PROJECT_ID/$SRR_ID/$mode"
    local lock="$CACHE/trimmed_fastq_manifest/$PROJECT_ID/${SRR_ID}_${mode}.lock"
    mkdir -p "$trim_dir" "$(dirname "$lock")"

    exec {trim_fd}>"$lock"
    flock "$trim_fd"

    local out1="$trim_dir/${SRR_ID}_${mode}_1.fastq.gz"
    local out2="$trim_dir/${SRR_ID}_${mode}_2.fastq.gz"
    if [ -s "$out1" ] && [ -s "$out2" ] && [ "$(stat -c%s "$out1")" -gt 1024 ] && [ "$(stat -c%s "$out2")" -gt 1024 ]; then
        append_row "$PREP_LOG" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$MODE	trim	cached	$TRIM_THREADS	0	0	$out1;$out2"
        printf "%s\t%s\n" "$out1" "$out2"
        flock -u "$trim_fd"
        return 0
    fi

    local trim_args
    if [ "$mode" = "adapter_only" ]; then
        trim_args="ILLUMINACLIP:$ADAPTER:2:30:10"
    else
        local pval=${mode#P}
        trim_args="ILLUMINACLIP:$ADAPTER:2:30:10 SLIDINGWINDOW:4:$pval LEADING:$pval TRAILING:$pval"
    fi

    local tmp_tag="${SLURM_JOB_ID:-manual}_${TASK_ID}_$$"
    local tmp1="$trim_dir/${SRR_ID}_${mode}_1.${tmp_tag}.tmp.fastq.gz"
    local tmp2="$trim_dir/${SRR_ID}_${mode}_2.${tmp_tag}.tmp.fastq.gz"
    local unp1="$trim_dir/${SRR_ID}_${mode}_1_unpaired.${tmp_tag}.tmp.fastq.gz"
    local unp2="$trim_dir/${SRR_ID}_${mode}_2_unpaired.${tmp_tag}.tmp.fastq.gz"
    local log="$ROOT/logs/${TASK_ID}_${SRR_ID}_${mode}_trimmomatic.log"
    local start end exit_code

    rm -f "$tmp1" "$tmp2" "$unp1" "$unp2"
    start=$(date +%s)
    set +e
    "$TRIMMOMATIC" -Xmx120g PE -threads "$TRIM_THREADS" -phred33 \
        "$in1" "$in2" "$tmp1" "$unp1" "$tmp2" "$unp2" \
        $trim_args > "$log" 2>&1
    exit_code=$?
    set -e
    end=$(date +%s)
    rm -f "$unp1" "$unp2"

    if [ "$exit_code" -ne 0 ] || [ ! -s "$tmp1" ] || [ ! -s "$tmp2" ] || [ "$(stat -c%s "$tmp1" 2>/dev/null || echo 0)" -le 1024 ] || [ "$(stat -c%s "$tmp2" 2>/dev/null || echo 0)" -le 1024 ]; then
        rm -f "$tmp1" "$tmp2"
        append_row "$PREP_LOG" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$MODE	trim	failed	$TRIM_THREADS	$((end - start))	$exit_code	$log"
        flock -u "$trim_fd"
        exit 1
    fi

    mv -f "$tmp1" "$out1"
    mv -f "$tmp2" "$out2"

    append_row "$PREP_LOG" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$MODE	trim	created	$TRIM_THREADS	$((end - start))	0	$log"
    printf "%s\t%s\n" "$out1" "$out2"
    flock -u "$trim_fd"
}

parse_star() {
    python3 - "$1" <<'PY'
import re
import sys
from pathlib import Path
text = Path(sys.argv[1]).read_text(errors="replace") if Path(sys.argv[1]).exists() else ""
def val(label):
    m = re.search(rf"{re.escape(label)}\s*\|\s*(.+)", text)
    return m.group(1).strip() if m else ""
print("\t".join([
    val("Uniquely mapped reads %"),
    val("Uniquely mapped reads number"),
    val("Number of reads mapped to multiple loci"),
    val("Number of input reads"),
]))
PY
}

parse_hisat2() {
    python3 - "$1" <<'PY'
import re
import sys
from pathlib import Path
text = Path(sys.argv[1]).read_text(errors="replace") if Path(sys.argv[1]).exists() else ""
rate = re.search(r"([\d.]+%)\s+overall alignment rate", text)
uniq = re.search(r"(\d+)\s+\([^)]+\)\s+aligned concordantly exactly 1 time", text)
multi = re.search(r"(\d+)\s+\([^)]+\)\s+aligned concordantly >1 times", text)
unmap = re.search(r"(\d+)\s+\([^)]+\)\s+aligned concordantly 0 times", text)
print("\t".join([
    rate.group(1) if rate else "",
    uniq.group(1) if uniq else "",
    multi.group(1) if multi else "",
    unmap.group(1) if unmap else "",
]))
PY
}

run_alignment() {
    local in1=$1
    local in2=$2
    local start end exit_code log parsed status read_cmd
    read_cmd=$(read_command_for "$in1")
    start=$(date +%s)
    if [ "$ALIGNER" = "STAR" ]; then
        local prefix="$WORK/star_"
        log="$ROOT/logs/${TASK_ID}_${SRR_ID}_${MODE}_STAR_Log.final.out"
        set +e
        "$STAR_BIN" --runMode alignReads \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$in1" "$in2" \
            --readFilesCommand "$read_cmd" \
            --runThreadN "$CORES" \
            --outSAMtype None \
            --outFileNamePrefix "$prefix" \
            --outTmpDir "$WORK/star_tmp" \
            > "$ROOT/logs/${TASK_ID}_${SRR_ID}_${MODE}_STAR.stdout" \
            2> "$ROOT/logs/${TASK_ID}_${SRR_ID}_${MODE}_STAR.stderr"
        exit_code=$?
        set -e
        [ -s "${prefix}Log.final.out" ] && cp "${prefix}Log.final.out" "$log"
        parsed=$(parse_star "$log" 2>/dev/null || printf "\t\t\t")
    elif [ "$ALIGNER" = "HISAT2" ]; then
        log="$ROOT/logs/${TASK_ID}_${SRR_ID}_${MODE}_HISAT2.log"
        if [ ! -s "${HISAT2_INDEX}.1.ht2" ] && [ ! -s "${HISAT2_INDEX}.1.ht2l" ]; then
            append_row "$RESULTS" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$CONDITION	$CELL_LINE	$MODE	$ALIGNER	$LAYOUT	$TOTAL_READS_M	$FASTQ_MB	$BOWTIE_TRIMMED_SEC	$BOWTIE_ALIGNMENT_RATE	$CORES	0	2	skipped_hisat2_index_missing					$log"
            return 0
        fi
        set +e
        "$HISAT2_BIN" -p "$CORES" -x "$HISAT2_INDEX" -1 "$in1" -2 "$in2" -S /dev/null 2> "$log"
        exit_code=$?
        set -e
        parsed=$(parse_hisat2 "$log" 2>/dev/null || printf "\t\t\t")
    else
        log="$ROOT/logs/${TASK_ID}_${SRR_ID}_${MODE}_BOWTIE2.log"
        if [ ! -s "${BOWTIE2_INDEX}.1.bt2" ]; then
            append_row "$RESULTS" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$CONDITION	$CELL_LINE	$MODE	$ALIGNER	$LAYOUT	$TOTAL_READS_M	$FASTQ_MB	$BOWTIE_TRIMMED_SEC	$BOWTIE_ALIGNMENT_RATE	$CORES	0	2	skipped_bowtie2_index_missing					$log"
            return 0
        fi
        if ! BOWTIE2_RESOLVED=$(find_bowtie2); then
            append_row "$RESULTS" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$CONDITION	$CELL_LINE	$MODE	$ALIGNER	$LAYOUT	$TOTAL_READS_M	$FASTQ_MB	$BOWTIE_TRIMMED_SEC	$BOWTIE_ALIGNMENT_RATE	$CORES	0	127	failed_bowtie2_not_found					$log"
            return 127
        fi
        set +e
        "$BOWTIE2_RESOLVED" -p "$CORES" -x "$BOWTIE2_INDEX" -1 "$in1" -2 "$in2" -S /dev/null 2> "$log"
        exit_code=$?
        set -e
        parsed=$(parse_hisat2 "$log" 2>/dev/null || printf "\t\t\t")
    fi
    end=$(date +%s)

    status="ok"
    [ "$exit_code" -eq 0 ] || status="failed"
    append_row "$RESULTS" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$CONDITION	$CELL_LINE	$MODE	$ALIGNER	$LAYOUT	$TOTAL_READS_M	$FASTQ_MB	$BOWTIE_TRIMMED_SEC	$BOWTIE_ALIGNMENT_RATE	$CORES	$((end - start))	$exit_code	$status	$parsed	$log"
}

echo "[$(date -Is)] task=$TASK_ID $PROJECT_ID $SRR_ID $MODE $ALIGNER"
if [ ! -s "$R1_PATH" ] || [ ! -s "$R2_PATH" ]; then
    append_row "$PREP_LOG" "$(date -Is)	$TASK_ID	$PROJECT_ID	$SRR_ID	$MODE	raw_fastq	missing	0	0	1	$R1_PATH;$R2_PATH"
    exit 1
fi

mode_pair=$(ensure_trimmed_fastqs "$MODE" "$R1_PATH" "$R2_PATH")
mode1=$(printf "%s" "$mode_pair" | cut -f1)
mode2=$(printf "%s" "$mode_pair" | cut -f2)

run_alignment "$mode1" "$mode2"
echo "[$(date -Is)] done task=$TASK_ID"
