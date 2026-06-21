#!/bin/bash
# Script Name: pipeline.sh
# Author: Fabian Braeuer
# Date: <2026-01>
# Description:
#   Starts sub-process for one SRR analysis.

BASE="$1"
PROJECT_ID="$2"
SRR_ID="$3"
PHRED_SCORES_STR="$4"
TRIMMERS_STR="$5"

source $BASE/scripts/src_pipeline/funtion_library/mv_slurm_logs.sh

pipelineLog_filePath="$BASE/results/$PROJECT_ID/$SRR_ID/pipeline.log"
srrData_dirPath="$BASE/raw_data/$PROJECT_ID/$SRR_ID"
LOG="$pipelineLog_filePath"
exec 2>> "$LOG"

out_err_dirPath="$BASE/results/$PROJECT_ID/$SRR_ID/_out_err"
mkdir -p "$out_err_dirPath"

job1_id=$(sbatch $BASE/scripts/src_pipeline/get_SRR_data_and_FastQc.sh $BASE $PROJECT_ID $SRR_ID | awk '{print $4}')
echo "submitted get_SRR_data_and_FastQc.sh ($job1_id) for $SRR_ID at $(date +"%Y-%m-%dT%H:%M:%S%z")" >&2
mv_slurm_logs "$BASE/scripts/run_get_SRR_data_and_FastQc_$job1_id" "$out_err_dirPath/run_get_SRR_data_and_FastQc_$job1_id" "$job1_id"

num_fastq_gz=$(find $srrData_dirPath -type f -name "*.fastq.gz" -size +0c | wc -l)

echo "num_fastq_gz value is: $num_fastq_gz" >&2

# job2_id=$(sbatch $BASE/scripts/src_pipeline/trimming.sh $BASE $PROJECT_ID $SRR_ID $num_fastq_gz "$PHRED_SCORES_STR" "$TRIMMERS_STR" | awk '{print $4}')
job2_id=$(sbatch --dependency=afterok:${job1_id} $BASE/scripts/src_pipeline/trimming.sh $BASE $PROJECT_ID $SRR_ID $num_fastq_gz "$PHRED_SCORES_STR" "$TRIMMERS_STR" | awk '{print $4}')
echo "submitted trimming.sh ($job2_id) for $SRR_ID at $(date +"%Y-%m-%dT%H:%M:%S%z")" >&2
mv_slurm_logs "$BASE/scripts/run_trimming_$job2_id" "$out_err_dirPath/run_trimming_$job2_id" "$job2_id"

# job3_id=$(sbatch $BASE/scripts/src_pipeline/alignment.sh $BASE $PROJECT_ID $SRR_ID $num_fastq_gz "$PHRED_SCORES_STR" "$TRIMMERS_STR" | awk '{print $4}')
job3_id=$(sbatch --dependency=afterok:${job2_id} $BASE/scripts/src_pipeline/alignment.sh $BASE $PROJECT_ID $SRR_ID $num_fastq_gz "$PHRED_SCORES_STR" "$TRIMMERS_STR" | awk '{print $4}')
echo "submitted alignment.sh ($job3_id) for $SRR_ID at $(date +"%Y-%m-%dT%H:%M:%S%z")" >&2
mv_slurm_logs "$BASE/scripts/run_alignment_$job3_id" "$out_err_dirPath/run_alignment_$job3_id" "$job3_id"

# job4_id=$(sbatch $BASE/scripts/src_pipeline/counting.sh $BASE $PROJECT_ID $SRR_ID $num_fastq_gz "$PHRED_SCORES_STR" "$TRIMMERS_STR" | awk '{print $4}')
job4_id=$(sbatch --dependency=afterok:${job3_id} $BASE/scripts/src_pipeline/counting.sh $BASE $PROJECT_ID $SRR_ID $num_fastq_gz "$PHRED_SCORES_STR" "$TRIMMERS_STR" | awk '{print $4}')
echo "submitted counting.sh ($job4_id) for $SRR_ID at $(date +"%Y-%m-%dT%H:%M:%S%z")" >&2
mv_slurm_logs "$BASE/scripts/run_counting_$job4_id" "$out_err_dirPath/run_counting_$job4_id" "$job4_id"
