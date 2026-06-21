#!/usr/bin/env python3
import os
import glob
import csv
import subprocess

METADATA_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/deseq2_metadata/sample_sheets"
CONCORDANCE_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance"
LOSO_SH = "/pc2/users/o/omiks001/scripts/analysis/loso_job.sh"
MAX_CORES = 128

def get_total_valid_samples(csv_path):
    # Matches the exact length of condition_map in 06_run_de_gsea.py
    total_valid = 0
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cond = row.get("Condition") or row.get("condition") or row.get("disease_state")
            if cond and str(cond).lower() not in ["none", "na", "null", "unknown"]:
                total_valid += 1
    return total_valid

def get_completed_samples(project_id):
    tsv_path = os.path.join(CONCORDANCE_DIR, f"{project_id}_concordance.tsv")
    if not os.path.exists(tsv_path):
        return 0
    
    completed = set()
    try:
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if "SRR_ID" in row:
                    completed.add(row["SRR_ID"])
    except Exception:
        pass
    return len(completed)

def submit_all():
    csv_files = glob.glob(os.path.join(METADATA_DIR, "*.csv"))
    print(f"Found {len(csv_files)} metadata files.")
    
    submitted_jobs = 0
    
    for csv_file in csv_files:
        project_id = os.path.basename(csv_file).replace(".csv", "")
        
        total_valid = get_total_valid_samples(csv_file)
        if total_valid < 4:
            continue # Needs at least 4 valid samples total (2 per group)
            
        completed = get_completed_samples(project_id)
        
        remaining = total_valid - completed
        if remaining <= 0:
            print(f"{project_id}: {total_valid} total, {completed} completed. Skipping.")
            continue
            
        # Cap max python workers at 20, but request 2x cores for 8GB RAM per worker
        python_workers = min(20, remaining + 2)
        slurm_cores = python_workers * 2
        
        print(f"{project_id}: {remaining} samples remaining ({total_valid} total). Requesting {slurm_cores} cores for {python_workers} workers.")
        
        cmd = [
            "sbatch",
            f"-c", str(slurm_cores),
            f"--job-name=loso_{project_id}",
            LOSO_SH,
            project_id,
            str(python_workers)
        ]
        
        try:
            res = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"  -> {res.stdout.strip()}")
            submitted_jobs += 1
        except subprocess.CalledProcessError as e:
            print(f"  -> Failed to submit {project_id}: {e.stderr}")
            
    print(f"Successfully submitted {submitted_jobs} jobs.")

if __name__ == "__main__":
    submit_all()
