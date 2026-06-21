#!/usr/bin/env python3

import os
import re
from pathlib import Path
from datetime import datetime
import csv

FB_RESULTS_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/fb_results")
FB_ORIGINAL_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/fb/omiks_project/results")

def parse_time(t_str):
    return datetime.strptime(t_str, "%Y-%m-%dT%H:%M:%S%z")

def main():
    out_files = []
    for base_dir in [FB_RESULTS_DIR, FB_ORIGINAL_DIR]:
        if base_dir.exists():
            out_files.extend(list(base_dir.rglob("_out_err/run_*counting*.out")))
            
    print(f"Found {len(out_files)} counting logs.")
    
    timings = []
    
    for f in out_files:
        try:
            text = f.read_text(encoding="utf-8", errors="replace")
        except:
            continue
            
        m_start = re.search(r"start(?:ed|et) .*?counting\.sh .*?at (\S+)", text)
        m_end = re.search(r"end(?:ed|et) .*?counting\.sh .*?at (\S+)", text)
        
        if m_start and m_end:
            start_str = m_start.group(1)
            end_str = m_end.group(1)
            start_dt = parse_time(start_str)
            end_dt = parse_time(end_str)
            duration = int((end_dt - start_dt).total_seconds())
            
            modes_processed = []
            if "counted .bam in" in text:
                for line in text.splitlines():
                    if "counted .bam in" in line:
                        if "/raw_data/" in line or "_untrimmed" in line:
                            modes_processed.append("untrimmed")
                        else:
                            m_mode = re.search(r"trimmomatic_([a-zA-Z0-9]+)\s+with", line)
                            if m_mode:
                                mode = m_mode.group(1)
                                if mode == "adapter":
                                    mode = "adapter_only"
                                modes_processed.append(mode)
            else:
                name = f.name
                if "trmd_counting" in name:
                    if "P35" in str(f): mode = "P35"
                    elif "P20" in str(f): mode = "P20"
                    elif "P10" in str(f): mode = "P10"
                    elif "P5" in str(f): mode = "P5"
                    elif "adapter" in str(f): mode = "adapter_only"
                    else: mode = "unknown"
                    if mode != "unknown":
                        modes_processed.append(mode)
                else:
                    modes_processed.append("untrimmed")
                    
            if not modes_processed:
                continue
                
            parts = f.parts
            srr_id = "unknown"
            project_id = "unknown"
            for i, p in enumerate(parts):
                if p.startswith("SRR"):
                    srr_id = p
                    project_id = parts[i-1]
                    break
                    
            if srr_id == "unknown":
                continue
                
            dur_per_mode = duration // len(modes_processed)
            
            for m in modes_processed:
                timings.append({
                    "project_id": project_id,
                    "srr_id": srr_id,
                    "mode": m,
                    "start_time": start_str,
                    "end_time": end_str,
                    "duration_sec": dur_per_mode
                })
                
    out_file = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_alignment/records/fb_counting_timing.tsv")
    with open(out_file, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["project_id", "srr_id", "mode", "start_time", "end_time", "duration_sec"], delimiter="\t")
        writer.writeheader()
        for t in timings:
            writer.writerow(t)
            
    print(f"Wrote {len(timings)} timing records to {out_file}.")

if __name__ == "__main__":
    main()
