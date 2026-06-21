#!/bin/bash
# v2_master_node.sh
# Usage: sbatch v2_master_node.sh <PROJECT_ID> <SRR_ID>

#SBATCH -J v2_pipeline
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=4G
#SBATCH -p normal

# We have requested 48 cores at 4GB each = 192GB.
# STAR is skipped, so memory is used for Bowtie2+Trimmomatic+featureCounts. We explicitly restrict the pipelines 
# to use fewer threads (e.g. 36) so they have massive spare memory overhead.
export ACTIVE_THREADS=36

PROJECT_ID=$1
SRR_ID=$2
LAYOUT=${3:-PAIRED}
FASTQ_PATHS=$4

if [ -z "$FASTQ_PATHS" ]; then
    echo "No FASTQ paths provided. Fetching data for $SRR_ID..."
    FETCH_TARGET="/scratch/hpc-prf-omiks/ja/temp/$PROJECT_ID/$SRR_ID"
    mkdir -p "$FETCH_TARGET"
    FASTQ_PATHS=$(bash "$HOME/scripts/fetch_data.sh" "$SRR_ID" "$FETCH_TARGET")
    if [ $? -ne 0 ] || [ -z "$FASTQ_PATHS" ]; then
        echo "ERROR: Fetch failed for $SRR_ID"
        exit 1
    fi
fi

if [ -z "$PROJECT_ID" ] || [ -z "$SRR_ID" ]; then
    echo "ERROR: Parameter missing. Expecting <PROJECT_ID> <SRR_ID> [FETCH_MODE]"
    exit 1
fi

DB_FILE="$HOME/srr_queue.db"
MAX_FETCH_FAILURES=10
BOWTIE_INDEX=${5:-/scratch/hpc-prf-omiks/fb/omiks_project/resources/indexes/GRCh38_index}
GTF_FILE=${6:-/scratch/hpc-prf-omiks/ja/omiks_project/resources/indexes/GRCh38_annotation.gtf}
DIR="$HOME/scripts"

# Boot the pipeline environment explicitly so batch jobs do not depend on an
# interactive submit shell having already activated omiks_pipeline.
CONDA_SH="/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/etc/profile.d/conda.sh"
PIPELINE_ENV_BIN="/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin"

if [ -f "$CONDA_SH" ]; then
    . "$CONDA_SH"
fi
export PATH="$PIPELINE_ENV_BIN:$PATH"
if command -v conda >/dev/null 2>&1; then
    conda activate omiks_pipeline >/dev/null 2>&1 || true
fi

for tool in trimmomatic bowtie2 samtools fastqc featureCounts; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: required tool '$tool' is not available on PATH"
        exit 1
    fi
done

# Prevent SRA tools from writing to $HOME
export NCBI_SETTINGS="/scratch/hpc-prf-omiks/ja/.ncbi/user-settings.mkfg"
export TMPDIR="/scratch/hpc-prf-omiks/ja/tmp"

# --- Permanent output directories ---
FASTQC_DIR="/scratch/hpc-prf-omiks/ja/flattened_fastqc_raw/$PROJECT_ID"
COUNTS_DIR="/scratch/hpc-prf-omiks/ja/flattened_counts/$PROJECT_ID"
TRIM_STATS_DIR="/scratch/hpc-prf-omiks/ja/flattened_trimmomatic_stats/$PROJECT_ID"
BOWTIE_STATS_DIR="/scratch/hpc-prf-omiks/ja/flattened_bowtie2_stats/$PROJECT_ID"
TRIM_STATS_FILE="$TRIM_STATS_DIR/${SRR_ID}_trimmomatic_stats.tsv"
BOWTIE_STATS_FILE="$BOWTIE_STATS_DIR/${SRR_ID}_bowtie2_stats.tsv"

# Use scratch-backed storage for all job intermediates.
WORK_DIR="$TMPDIR/v2_job_${SLURM_JOB_ID}_${SRR_ID}"
echo "Using scratch work directory for temp storage: $WORK_DIR"

mkdir -p "$WORK_DIR"
mkdir -p "$TRIM_STATS_DIR" "$BOWTIE_STATS_DIR"

# --- Timing setup ---
TIMING_DIR="/scratch/hpc-prf-omiks/ja/flattened_timings/$PROJECT_ID"
mkdir -p "$TIMING_DIR"
TIMING_FILE="$TIMING_DIR/${SRR_ID}_timings.tsv"
PIPELINE_START=$SECONDS
echo -e "stage\tstart_epoch\tend_epoch\tduration_sec" > "$TIMING_FILE"
[ -f "$TRIM_STATS_FILE" ] || echo -e "project_id\tsrr_id\tmode\tlayout\tinput_reads\tsurviving\tsurviving_pct\tdropped\tdropped_pct" > "$TRIM_STATS_FILE"
[ -f "$BOWTIE_STATS_FILE" ] || echo -e "project_id\tsrr_id\tmode\tlayout\ttotal_reads\tpaired_reads\tunpaired_reads\taligned_exactly_1\taligned_gt1\taligned_0\tconcordant_exactly_1\tconcordant_gt1\tdiscordant_1\tpairs_0_concordant\toverall_alignment_rate" > "$BOWTIE_STATS_FILE"

timer_start() { eval "T_${1}=$SECONDS"; }
timer_end() {
  local name=$1
  local start_var="T_${name}"
  local start=${!start_var}
  local dur=$(( SECONDS - start ))
  echo -e "${name}\t${start}\t${SECONDS}\t${dur}" >> "$TIMING_FILE"
  echo "[TIMING] $name: ${dur}s"
}

# --- Checkpoint helper ---
stage_done() {
    # Check if any file matching the pattern exists
    ls $1 &>/dev/null
}

echo "=========================================================="
echo "V2 Pipeline with Checkpoint/Resume (STAR SKIPPED)"
echo "Hostname: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Target SRR: $SRR_ID ($PROJECT_ID)"
echo "Start Time: $(date)"
echo "=========================================================="

function trigger_rollback() {
    echo "PIPELINE FAILED: Triggering secure SQLite DB rollback."
    rm -rf "$WORK_DIR"
    sqlite3 "$DB_FILE" "UPDATE srr_queue
        SET pipeline_failures=COALESCE(pipeline_failures,0)+1,
            run_status='FETCHED',
            status='todo',
            slurm_job_id=NULL,
            last_updated=CURRENT_TIMESTAMP
        WHERE srr_id='$SRR_ID';"
    exit 1
}

function trigger_refetch() {
    local reason="${1:-invalid_input_files}"
    echo "INPUT INVALID: forcing refetch for $SRR_ID ($reason)."
    rm -rf "$WORK_DIR"
    rm -f "$FASTQC_DIR"/${SRR_ID}_*_fastqc.html "$FASTQC_DIR"/${SRR_ID}_*_fastqc.zip 2>/dev/null
    rm -rf "$FASTQC_DIR"/${SRR_ID}_*_fastqc 2>/dev/null
    sqlite3 "$DB_FILE" "UPDATE srr_queue
        SET failure_count=COALESCE(failure_count,0)+1,
            run_status=CASE
                WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'FETCH_ABANDONED'
                ELSE 'FETCH_FAILED'
            END,
            status=CASE
                WHEN COALESCE(failure_count,0)+1 >= $MAX_FETCH_FAILURES THEN 'abandoned'
                ELSE 'todo'
            END,
            local_fastq_path=NULL,
            slurm_job_id=NULL,
            last_updated=CURRENT_TIMESTAMP
        WHERE srr_id='$SRR_ID';"
    exit 1
}

# ===== STAGE 1: SRA CONVERSION (ON LUSTRE) =====
# SRA conversion extracts massive uncompressed FASTQ files, so keep it on scratch.
timer_start sra_conversion
echo -e "\n[$(date)] >>> STAGE 1: SRA CONVERSION (Lustre) <<<"

SRA_FILE=$(echo "$FASTQ_PATHS" | tr ' ' '\n' | grep '\.sra$' | head -1)
if [ -n "$SRA_FILE" ] && [ -f "$SRA_FILE" ]; then
    sqlite3 "$DB_FILE" "UPDATE srr_queue SET run_status='SRA_CONVERSION' WHERE srr_id='$SRR_ID';"
    SRA_DIR=$(dirname "$SRA_FILE")
    
    echo "Extracting SRA safely on Lustre ($SRA_DIR)..."
    fasterq-dump -e "$ACTIVE_THREADS" --split-files --outdir "$SRA_DIR" "$SRA_FILE"
    FQD_EXIT=$?

    if [ $FQD_EXIT -ne 0 ]; then
        echo "fasterq-dump failed (exit $FQD_EXIT), falling back to fastq-dump..."
        fastq-dump --split-files "$SRA_FILE" -O "$SRA_DIR/"
        if [ $? -ne 0 ]; then
            echo "ERROR: Both fasterq-dump and fastq-dump failed for $SRA_FILE"
            rm -f "$SRA_FILE"
            trigger_refetch "sra_conversion_failed"
        fi
    fi

    # Immediately delete the SRA file once converted to save space
    rm -f "$SRA_FILE"

    for f in "$SRA_DIR"/*.fastq; do
        if [ -f "$f" ]; then
            if command -v pigz &>/dev/null; then
                pigz -p "$ACTIVE_THREADS" "$f"
            else
                gzip "$f"
            fi
        fi
    done

    FASTQ_PATHS=$(ls -1d "$SRA_DIR"/*.fastq.gz 2>/dev/null | tr '\n' ' ' | xargs)
    echo "Converted to: $FASTQ_PATHS"

    if [ -z "$FASTQ_PATHS" ]; then
        echo "ERROR: No .fastq.gz files after SRA conversion"
        trigger_refetch "sra_conversion_produced_no_fastq"
    fi

    FQ_BYTES=$(du -cb $FASTQ_PATHS | awk '/total/ {print $1}')
    sqlite3 "$DB_FILE" "UPDATE srr_queue SET fastq_size_bytes=$FQ_BYTES WHERE srr_id='$SRR_ID';"
fi
timer_end sra_conversion

# ===== STAGE 1b: LOCAL COPY (SKIPPED) =====
# Base FASTQ.gz files stay on scratch and stream into each mode.

# ===== STAGE 1c: FASTQ INTEGRITY CHECK =====
echo -e "\n[$(date)] >>> STAGE 1c: FASTQ INTEGRITY CHECK <<<"
ORIG_FASTQ_PATHS="$4"  # original paths passed to script (on /scratch)
for fq in $FASTQ_PATHS; do
    if [ ! -s "$fq" ]; then
        echo "ERROR: $fq is empty (0 bytes). Deleting source files."
        for src in $ORIG_FASTQ_PATHS; do rm -f "$src"; done
        trigger_refetch "fastq_empty"
    fi
    gzip -t "$fq" 2>/dev/null
    if [ $? -ne 0 ]; then
        echo "ERROR: $fq is corrupt (gzip -t failed). Deleting source files."
        for src in $ORIG_FASTQ_PATHS; do rm -f "$src"; done
        trigger_refetch "fastq_gzip_corrupt"
    fi
done
echo "All FASTQ files passed integrity check."

# ===== STAGE 1d: RECORD FASTQ METADATA FOR THROUGHPUT ANALYSIS =====
FASTQ_BYTES=0
for fq in $FASTQ_PATHS; do
    FQ_SIZE=$(stat --printf='%s' "$fq" 2>/dev/null || stat -f '%z' "$fq" 2>/dev/null || echo 0)
    FASTQ_BYTES=$((FASTQ_BYTES + FQ_SIZE))
done
FASTQ_MB=$((FASTQ_BYTES / 1048576))
echo -e "fastq_size_bytes\t$FASTQ_BYTES" >> "$TIMING_FILE"
echo -e "fastq_size_mb\t$FASTQ_MB" >> "$TIMING_FILE"
echo "FASTQ total size: ${FASTQ_MB} MB (compressed)"

# ===== STAGE 2: FASTQC =====
mkdir -p "$FASTQC_DIR"
if stage_done "$FASTQC_DIR/${SRR_ID}_*_fastqc.html"; then
    echo -e "\n[SKIP] FastQC already done for $SRR_ID"
else
    timer_start fastqc
    echo -e "\n[$(date)] >>> STAGE 2: RAW FASTQC <<<"
    sqlite3 "$DB_FILE" "UPDATE srr_queue SET run_status='FASTQC (RAW)' WHERE srr_id='$SRR_ID';"

    FQC_PIDS=()
    for fq in $FASTQ_PATHS; do
        echo "Running FastQC on: $fq"
        fastqc -quiet -t 8 "$fq" -o "$FASTQC_DIR" --extract &
        FQC_PIDS+=($!)
    done
    for pid in "${FQC_PIDS[@]}"; do
        wait $pid
        if [ $? -ne 0 ]; then trigger_rollback; fi
    done
    echo "FastQC exported to $FASTQC_DIR"
    timer_end fastqc
fi

# Extract total reads + read length from FastQC output and record in timing file
for fqc_dir in "$FASTQC_DIR"/${SRR_ID}*_fastqc/; do
    [ -d "$fqc_dir" ] || continue
    DATA_FILE="$fqc_dir/fastqc_data.txt"
    [ -f "$DATA_FILE" ] || continue
    TOTAL_SEQS=$(awk -F'\t' '/^Total Sequences/ {print $2; exit}' "$DATA_FILE")
    SEQ_LEN=$(awk -F'\t' '/^Sequence length/ {print $2; exit}' "$DATA_FILE")
    if [ -n "$TOTAL_SEQS" ]; then
        echo -e "total_reads\t$TOTAL_SEQS" >> "$TIMING_FILE"
        echo -e "read_length\t$SEQ_LEN" >> "$TIMING_FILE"
        echo "Recorded: total_reads=$TOTAL_SEQS, read_length=$SEQ_LEN"
        break  # Only need one read file's stats
    fi
done

# ===== STAGE 3: METADATA DETECTION =====
echo -e "\n[$(date)] >>> STAGE 3: FASTQ METADATA DETECTION <<<"
METADATA=$(bash "$DIR/check_fastq_metadata.sh" $FASTQ_PATHS)
if [ $? -ne 0 ]; then trigger_rollback; fi
eval "$METADATA"
echo "Parsed Layout: $LAYOUT, Encoding: $PHRED"

# ===== STAGE 4: EVAL MODES (trim/align/count × 6) =====
MODES=("untrimmed" "adapter_only" "P5" "P10" "P20" "P35")
mkdir -p "$COUNTS_DIR"

sqlite3 "$DB_FILE" "UPDATE srr_queue SET run_status='EVALUTING ALL MODES' WHERE srr_id='$SRR_ID';"

timer_start eval_modes
SKIPPED_MODES=0
for MODE_RAW in "${MODES[@]}"; do
    # Determine the count output filename for this mode
    if [ "$MODE_RAW" = "untrimmed" ]; then
        MODE_NAME="untrmd_${SRR_ID}"
    elif [ "$MODE_RAW" = "adapter_only" ]; then
        MODE_NAME="${SRR_ID}_trimmomatic_adapter"
    else
        MODE_NAME="${SRR_ID}_trimmomatic_$MODE_RAW"
    fi

    COUNT_FILE="$COUNTS_DIR/${MODE_NAME}_fC.txt.gz"

    if [ -f "$COUNT_FILE" ]; then
        echo "[SKIP] Mode $MODE_RAW already done ($COUNT_FILE exists)"
        SKIPPED_MODES=$((SKIPPED_MODES + 1))
        continue
    fi

    echo -e "\n[$(date)] >>> STAGE 4: MODE $MODE_RAW [Threads: $ACTIVE_THREADS (Slurm: $SLURM_CPUS_PER_TASK)] <<<"
    MODE_START=$SECONDS

    # Base FASTQs stream from Lustre while Trimmomatic, Bowtie2, samtools sort,
    # and featureCounts intermediates live in scratch for this mode.
    MODE_WORK_DIR="$WORK_DIR/$MODE_RAW"
    mkdir -p "$MODE_WORK_DIR"

    # TRIMMING
    T_TRIM=$SECONDS
    if [ "$MODE_RAW" = "untrimmed" ]; then
        TRIM_PATHS="$FASTQ_PATHS"
    else
        export TRIM_STATS_FILE PROJECT_ID
        TRIM_PATHS=$(bash "$DIR/v2_trim.sh" "$MODE_RAW" "$SRR_ID" "$LAYOUT" "$PHRED" "$MODE_WORK_DIR" $FASTQ_PATHS)
        if [ $? -ne 0 ]; then trigger_rollback; fi
    fi
    echo -e "trim_${MODE_RAW}\t${T_TRIM}\t${SECONDS}\t$(( SECONDS - T_TRIM ))" >> "$TIMING_FILE"

    # ALIGNMENT
    T_ALIGN=$SECONDS
    export BOWTIE_STATS_FILE PROJECT_ID MODE_NAME
    BAM_PATH=$(bash "$DIR/v2_align.sh" "$SRR_ID" "$LAYOUT" "$PHRED" "$BOWTIE_INDEX" "$MODE_WORK_DIR" $TRIM_PATHS)
    if [ $? -ne 0 ]; then trigger_rollback; fi
    echo -e "align_${MODE_RAW}\t${T_ALIGN}\t${SECONDS}\t$(( SECONDS - T_ALIGN ))" >> "$TIMING_FILE"

    # COUNTING (writes to permanent storage)
    T_COUNT=$SECONDS
    bash "$DIR/v2_count.sh" "$MODE_NAME" "$PROJECT_ID" "$BAM_PATH" "$GTF_FILE" "$LAYOUT"
    if [ $? -ne 0 ]; then trigger_rollback; fi
    echo -e "count_${MODE_RAW}\t${T_COUNT}\t${SECONDS}\t$(( SECONDS - T_COUNT ))" >> "$TIMING_FILE"

    echo -e "mode_${MODE_RAW}_total\t${MODE_START}\t${SECONDS}\t$(( SECONDS - MODE_START ))" >> "$TIMING_FILE"

    # CLEANUP mode temp
    rm -rf "$MODE_WORK_DIR"
done

if [ "$SKIPPED_MODES" -eq 6 ]; then
    echo "[$(date)] All 6 modes already completed"
else
    echo "[$(date)] ALL MODES COMPLETED ($SKIPPED_MODES skipped, $((6 - SKIPPED_MODES)) ran)"
fi
timer_end eval_modes

# ===== STAR + DOWNSTREAM SKIPPED =====
echo -e "\n[$(date)] Skipping STAR mapping, sex calling, and virus profiling (not needed for current run)."

# ===== DONE =====
PIPELINE_DUR=$(( SECONDS - PIPELINE_START ))
echo -e "pipeline_total\t${PIPELINE_START}\t${SECONDS}\t${PIPELINE_DUR}" >> "$TIMING_FILE"
echo -e "\n[$(date)] PIPELINE COMPLETED SUCCESSFULLY (total: ${PIPELINE_DUR}s)"
echo "Timings written to: $TIMING_FILE"
rm -rf "$WORK_DIR"
sqlite3 "$DB_FILE" "UPDATE srr_queue SET run_status='DONE_NO_STAR', status='done', last_updated=CURRENT_TIMESTAMP WHERE srr_id='$SRR_ID';"
exit 0
