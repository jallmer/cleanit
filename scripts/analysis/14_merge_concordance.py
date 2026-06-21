#!/usr/bin/env python3
import os
import glob
import pandas as pd

CONCORDANCE_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance"
OUT_FILE = os.path.join(CONCORDANCE_DIR, "all_concordance.tsv")

def merge_concordance():
    # Find all individual project TSV files
    files = glob.glob(os.path.join(CONCORDANCE_DIR, "*_concordance.tsv"))
    
    # Exclude the output file itself and bio_concordance (legacy)
    files = [f for f in files if not f.endswith("all_concordance.tsv") and not f.endswith("bio_concordance.tsv")]
    
    print(f"Found {len(files)} individual concordance files.")
    
    all_dfs = []
    for fp in files:
        try:
            df = pd.read_csv(fp, sep="\t")
            all_dfs.append(df)
        except Exception as e:
            print(f"Skipping {fp}: {e}")
            
    if not all_dfs:
        print("No data found to merge.")
        return
        
    master_df = pd.concat(all_dfs, ignore_index=True)
    
    # Filter out rows where NA is present (the failed edge cases)
    # We drop rows where U_rho_gene is NA
    initial_len = len(master_df)
    master_df = master_df.dropna(subset=["U_rho_gene"])
    
    # Drop any duplicates that occurred due to overlapping job submissions
    master_df = master_df.drop_duplicates(subset=["SRR_ID", "project_id"])
    
    print(f"Merged {initial_len} total rows. Filtered down to {len(master_df)} valid unique completed samples.")
    
    master_df.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"Saved master table to {OUT_FILE}")

if __name__ == "__main__":
    merge_concordance()
