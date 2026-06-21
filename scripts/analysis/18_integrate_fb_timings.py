#!/usr/bin/env python3

import os
import csv
from pathlib import Path
import pandas as pd

BASE_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks")
JA_TIMINGS_DIR = BASE_DIR / "ja/flattened_timings"

FB_TRIM_TIMING = BASE_DIR / "ja/final_trimming/records/fb_timing_table.tsv"
FB_ALIGN_TIMING = BASE_DIR / "ja/final_alignment/records/fb_alignment_timing.tsv"
FB_COUNT_TIMING = BASE_DIR / "ja/final_alignment/records/fb_counting_timing.tsv"
COMBINED_STATS = BASE_DIR / "ja/final_alignment/records/combined_alignment_stats.tsv"

def load_fb_timings(filepath):
    if not filepath.exists():
        return {}
    df = pd.read_csv(filepath, sep="\t")
    if 'end_time' in df.columns:
        df['end_time_dt'] = pd.to_datetime(df['end_time'], format="%Y-%m-%dT%H:%M:%S%z", errors='coerce')
        df = df.sort_values(by="end_time_dt", ascending=False)
        df = df.drop_duplicates(subset=["project_id", "srr_id", "mode"], keep="first")
    
    # Store clean modes
    result = {}
    for _, row in df.iterrows():
        srr = row["srr_id"]
        mode = row["mode"]
        clean_mode = mode
        if "untrimmed" in mode or "untrmd" in mode:
            clean_mode = "untrimmed"
        elif "adapter" in mode:
            clean_mode = "adapter_only"
        elif "P5" in mode: clean_mode = "P5"
        elif "P10" in mode: clean_mode = "P10"
        elif "P20" in mode: clean_mode = "P20"
        elif "P35" in mode: clean_mode = "P35"
        
        result[(srr, clean_mode)] = row["duration_sec"]
    return result

def load_ja_timings(filepath):
    timings = {}
    if filepath.exists():
        with open(filepath, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                stage = row.get("stage", "")
                dur = row.get("duration_sec", "")
                timings[stage] = int(float(dur)) if dur and dur != "NA" else 0
    return timings

def save_ja_timings(filepath, timings):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["stage", "start_epoch", "end_epoch", "duration_sec"])
        for k, v in timings.items():
            writer.writerow([k, 0, 0, v])

def main():
    print("Loading datasets...")
    trim_dict = load_fb_timings(FB_TRIM_TIMING)
    align_dict = load_fb_timings(FB_ALIGN_TIMING)
    count_dict = load_fb_timings(FB_COUNT_TIMING)
    
    combined_df = pd.read_csv(COMBINED_STATS, sep="\t")
    fb_df = combined_df[combined_df["source"] == "fb"]
    
    srr_modes = {}
    for _, row in fb_df.iterrows():
        srr = row["srr_id"]
        proj = row["project_id"]
        mode = row["mode"]
        
        clean_mode = mode
        if "untrimmed" in mode or "untrmd" in mode:
            clean_mode = "untrimmed"
        elif "adapter" in mode:
            clean_mode = "adapter_only"
        elif "P5" in mode: clean_mode = "P5"
        elif "P10" in mode: clean_mode = "P10"
        elif "P20" in mode: clean_mode = "P20"
        elif "P35" in mode: clean_mode = "P35"
        
        if srr not in srr_modes:
            srr_modes[srr] = {"project": proj, "modes": []}
        srr_modes[srr]["modes"].append((mode, clean_mode))
        
    print(f"Integrating timings for {len(srr_modes)} SRRs...")
    
    updates = 0
    created = 0
    
    for srr, info in srr_modes.items():
        proj = info["project"]
        timing_file = JA_TIMINGS_DIR / proj / f"{srr}_timings.tsv"
        
        exists = timing_file.exists()
        timings = load_ja_timings(timing_file)
        
        for raw_mode, clean_mode in info["modes"]:
            t_trim = trim_dict.get((srr, clean_mode), 0)
            t_align = align_dict.get((srr, clean_mode), 0)
            t_count = count_dict.get((srr, clean_mode), 0)
            
            if clean_mode == "untrimmed":
                t_trim = 0  
                
            total = t_trim + t_align + t_count
            
            timings[f"trim_{clean_mode}"] = t_trim
            timings[f"align_{clean_mode}"] = t_align
            timings[f"count_{clean_mode}"] = t_count
            timings[f"mode_{clean_mode}_total"] = total
            
        if "pipeline_total" not in timings:
            timings["pipeline_total"] = sum([v for k, v in timings.items() if k.startswith("mode_") and k.endswith("_total")]) + timings.get("fastqc", 0) + timings.get("sra_conversion", 0)
            
        save_ja_timings(timing_file, timings)
        
        if exists:
            updates += 1
        else:
            created += 1
            
    print(f"Integration complete.")
    print(f"  Updated existing SRR files: {updates}")
    print(f"  Created new SRR files: {created}")

if __name__ == "__main__":
    main()
