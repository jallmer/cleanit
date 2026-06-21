#!/bin/bash
set -e

echo "Starting QC fitting pipeline..."

echo "1. Merging Concordance Results..."
/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/python3 /pc2/users/o/omiks001/scripts/analysis/14_merge_concordance.py

echo "2. Classifying Optimal Trimming (t*)..."
/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/python3 /pc2/users/o/omiks001/scripts/analysis/07_classify_trimming.py

echo "3. Extracting FastQC stats..."
/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/python3 /pc2/users/o/omiks001/scripts/analysis/01_extract_fastqc_stats.py

echo "4. Building Sample Feature Table..."
/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/python3 /pc2/users/o/omiks001/scripts/analysis/13_build_sample_feature_table.py

echo "5. Fitting QC Models..."
/pc2/users/o/omiks001/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/python3 /pc2/users/o/omiks001/scripts/analysis/08_fit_qc_model.py

echo "QC fitting pipeline complete!"

