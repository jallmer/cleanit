#!/bin/bash
# Script Name: start_missing_pipeline.sh
# Description: Custom runner that scans the flattened count directory and only executes the 
# pipeline on SRRs that are genuinely missing output text files.

# Ensure we use your personal JA working directory
BASE="/scratch/hpc-prf-omiks/ja/omiks_project"
JA_COUNTS="/scratch/hpc-prf-omiks/ja/flattened_counts"

PHRED_SCORES_STR="35"
TRIMMERS_STR="trimmomatic"

data_dirPath="$BASE/raw_data"
results_dirPath="$BASE/results"

ERROR_LOG="$BASE/scripts/start_pipeline_error.log"
: > "$ERROR_LOG"
exec 2>> "$ERROR_LOG"

project_IDs_filtered_filePath="$BASE/resources/project_IDs_filtered.tsv"        

if [[ ! -s "$project_IDs_filtered_filePath" ]]; then
    echo "ERROR: $project_IDs_filtered_filePath not found. Please pull the repository into your /ja/ workspace." >&2
    exit 1
fi

mkdir -p "$data_dirPath" "$results_dirPath"

project_counter=0
while IFS=$'\t' read -r project_ID srr_count_anticipated; do
    [[ "$project_ID" =~ ^# ]] && continue
    
    ((project_counter++))

    projectData_dirPath="$data_dirPath/$project_ID"
    projectResults_dirPath="$results_dirPath/$project_ID"
    mkdir -p "$projectData_dirPath" "$projectResults_dirPath"
    
    srr_IDs_filePath="$BASE/resources/bioproject_SRR_IDs/${project_ID}_SRR_IDs.txt"
    if [ ! -f "$srr_IDs_filePath" ]; then continue; fi

    mapfile -t SRR_IDs_arr < "$srr_IDs_filePath"
    srr_counter=0
    for srr_ID in "${SRR_IDs_arr[@]}"; do
        srr_ID="${srr_ID//$'\r'/}" 

        # --- MISSING CHECK LOGIC ---
        # 6 Expected files per SRR
        missing=0
        for suffix in "untrmd_${srr_ID}" "${srr_ID}_trimmomatic_adapter" "${srr_ID}_trimmomatic_P5" "${srr_ID}_trimmomatic_P10" "${srr_ID}_trimmomatic_P20" "${srr_ID}_trimmomatic_P35"; do
            if [ ! -f "$JA_COUNTS/$project_ID/${suffix}_fC.txt.gz" ]; then
                missing=1
                break
            fi
        done
        
        if [ $missing -eq 0 ]; then
            # echo "Skipping $srr_ID (fully processed and flattened)"
            ((srr_counter++))
            continue
        fi
        
        echo "Submitting missing job for $srr_ID..."

        mkdir -p "$projectData_dirPath/$srr_ID" "$projectResults_dirPath/$srr_ID" 

        bash "$BASE/scripts/src_pipeline/pipeline.sh" \
            "$BASE" \
            "$project_ID" \
            "$srr_ID" \
            "$PHRED_SCORES_STR" \
            "$TRIMMERS_STR"
            
        status_pipeline=$?
        if (( status_pipeline != 0 )); then
            echo "ERROR: pipeline.sh failed to start for SRR $srr_ID" >&2
            break
        fi
        ((srr_counter++))
        sleep 0.1
    done

    (( project_counter == 1 )) && printf "%-15s\t%-20s\t%-20s\t%-20s\n" "project ID" "anticipated SRRs" "processed SRRs"
    mismatch=$([ "$srr_count_anticipated" -ne "$srr_counter" ] && echo "<- MISMATCH" || echo "")
    printf "%-15s\t%-20d\t%-20d\t%-20s\n" "$project_ID" "$srr_count_anticipated" "$srr_counter" "$mismatch"

done < "$project_IDs_filtered_filePath"

echo "Missing project evaluation completed."
