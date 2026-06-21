#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import gzip
import shutil
import sys

def compress_file(src_path, dest_path):
    if dest_path.exists():
        return False
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(src_path, 'rb') as f_in:
        with gzip.open(dest_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return True

def main():
    base_dir = Path("/pc2/users/o/omiks001")
    ja_final_counts = base_dir / "hpc-prf-omiks/ja/flattened_counts"
    fb_results_1 = base_dir / "hpc-prf-omiks/fb/omiks_project/results"
    fb_results_2 = base_dir / "hpc-prf-omiks/ja/final_trimming/records/fb_results"
    
    combined_stats_path = base_dir / "hpc-prf-omiks/ja/final_alignment/records/combined_alignment_stats.tsv"
    
    if not combined_stats_path.exists():
        print("Combined stats not found.")
        sys.exit(1)
        
    df = pd.read_csv(combined_stats_path, sep="\t")
    fb_df = df[df['source'] == 'fb']
    
    print(f"Total fb records to process: {len(fb_df)}")
    
    success_count = 0
    missing_count = 0
    skip_count = 0
    
    for _, row in fb_df.iterrows():
        project_id = row['project_id']
        srr_id = row['srr_id']
        mode = row['mode']
        
        if mode == "untrimmed":
            prefix = f"untrmd_{srr_id}_fC.txt"
        else:
            if mode == "adapter_only":
                fb_mode_str = "adapter"
            else:
                fb_mode_str = mode
                
            if fb_mode_str.startswith("SRR") and "_trimmomatic_" in fb_mode_str:
                prefix = f"{fb_mode_str}_fC.txt"
            else:
                prefix = f"{srr_id}_trimmomatic_{fb_mode_str}_fC.txt"
                
        for suffix in ["", ".summary"]:
            filename = prefix + suffix
            target_filename = filename + ".gz"
            target_path = ja_final_counts / project_id / target_filename
            
            found_src = None
            for results_dir in [fb_results_1, fb_results_2]:
                srr_dir = results_dir / project_id / srr_id
                if not srr_dir.exists():
                    continue
                for candidate in srr_dir.rglob(filename):
                    found_src = candidate
                    break
                if found_src:
                    break
                    
            if found_src:
                if target_path.exists():
                    skip_count += 1
                else:
                    if compress_file(found_src, target_path):
                        success_count += 1
            else:
                missing_count += 1

    print(f"Successfully migrated: {success_count} files")
    print(f"Skipped (already existed): {skip_count} files")
    print(f"Missing files: {missing_count} files")

if __name__ == "__main__":
    main()
