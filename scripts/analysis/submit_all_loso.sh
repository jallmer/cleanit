#!/bin/bash
METADATA_DIR="/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/deseq2_metadata/sample_sheets"
LOGS_DIR="/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/logs"

mkdir -p "$LOGS_DIR"

echo "Submitting LOSO jobs for all suitable bioprojects..."
for csv in "$METADATA_DIR"/*.csv; do
    if [ -f "$csv" ]; then
        project_id=$(basename "$csv" .csv)
        sbatch /pc2/users/o/omiks001/scripts/analysis/loso_job.sh "$project_id"
    fi
done
echo "All jobs submitted!"
