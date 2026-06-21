#!/usr/bin/env python3
"""
Merge student fb_results trimmomatic stats and extract timings.

Reads logs from:
  /pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/fb_results/<PROJECT>/<SRR>/_out_err/
And:
  /pc2/users/o/omiks001/hpc-prf-omiks/fb/omiks_project/results/<PROJECT>/<SRR>/_out_err/

Writes to:
  /pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/flattened_trimmomatic_stats/<PROJECT>/<SRR>_trimmomatic_stats.tsv
  /pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/fb_timing_table.tsv
"""

import sqlite3
import os
import glob
import re
import csv
from datetime import datetime
from pathlib import Path

DB_FILE = os.path.expanduser("~/srr_queue.db")
FB_RESULTS_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/fb_results")
FB_ORIGINAL_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/fb/omiks_project/results")
TRIM_BASE = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/flattened_trimmomatic_stats")
TIMING_FILE = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/fb_timing_table.tsv")

TRIM_HEADER = [
    "project_id",
    "srr_id",
    "mode",
    "layout",
    "input_reads",
    "surviving",
    "surviving_pct",
    "dropped",
    "dropped_pct",
]

TIMING_HEADER = ["project_id", "srr_id", "mode", "start_time", "end_time", "duration_sec"]

def get_srrs():
    conn = sqlite3.connect(DB_FILE)
    cur = conn.cursor()
    cur.execute("SELECT project_id, srr_id FROM srr_queue")
    srrs = cur.fetchall()
    conn.close()
    return srrs

def parse_trimmomatic_stats(err_text, project_id, srr_id):
    # Try to find mode from command line inside err_text
    # -threads 8 ... SRR12676593_trimmomatic_P35_trmd.fastq.gz
    m = re.search(f"{srr_id}_trimmomatic_([a-zA-Z0-9]+)_trmd", err_text)
    mode = m.group(1) if m else "unknown"
    if mode == "unknown":
        # Check if mode is 'adapter' without trmd in some cases
        if f"{srr_id}_trimmomatic_adapter" in err_text:
            mode = "adapter_only" # To match flattened format
        else:
            # Fallback
            for md in ["P35", "P20", "P10", "P5", "adapter"]:
                if f"_{md}_" in err_text or f"_{md}." in err_text:
                    mode = "adapter_only" if md == "adapter" else md
                    break
    if mode == "adapter":
        mode = "adapter_only"

    # Single-end
    se_re = re.search(r"Input Reads:\s+(\d+)\s+Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Dropped:\s+(\d+)\s+\(([\d.]+)%\)", err_text)
    if se_re:
        return {
            "project_id": project_id,
            "srr_id": srr_id,
            "mode": mode,
            "layout": "SE",
            "input_reads": se_re.group(1),
            "surviving": se_re.group(2),
            "surviving_pct": se_re.group(3),
            "dropped": se_re.group(4),
            "dropped_pct": se_re.group(5),
        }
    
    # Paired-end
    pe_re = re.search(r"Input Read Pairs:\s+(\d+)\s+Both Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Forward Only Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Reverse Only Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Dropped:\s+(\d+)\s+\(([\d.]+)%\)", err_text)
    if pe_re:
        input_reads = int(pe_re.group(1))
        both = int(pe_re.group(2))
        forward = int(pe_re.group(4))
        reverse = int(pe_re.group(6))
        dropped = int(pe_re.group(8))
        surviving = both + forward + reverse
        surviving_pct = f"{(surviving / input_reads * 100) if input_reads else 0:.2f}"
        return {
            "project_id": project_id,
            "srr_id": srr_id,
            "mode": mode,
            "layout": "PE",
            "input_reads": str(input_reads),
            "surviving": str(surviving),
            "surviving_pct": surviving_pct,
            "dropped": str(dropped),
            "dropped_pct": pe_re.group(9),
        }
    
    return None

def parse_timing(out_text, mode):
    # startet trimmomatic.sh for SRR12676593 at 2026-02-24T15:18:26+0100
    # ended trimmomatic.sh for SRR12676593 at 2026-02-24T15:19:54+0100
    start_m = re.search(r"startet\s+.*?at\s+(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}[+-]\d{4})", out_text)
    end_m = re.search(r"ended\s+.*?at\s+(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}[+-]\d{4})", out_text)
    
    start_time = start_m.group(1) if start_m else ""
    end_time = end_m.group(1) if end_m else ""
    duration_sec = ""

    if start_time and end_time:
        try:
            fmt = "%Y-%m-%dT%H:%M:%S%z"
            s_dt = datetime.strptime(start_time, fmt)
            e_dt = datetime.strptime(end_time, fmt)
            duration_sec = str(int((e_dt - s_dt).total_seconds()))
        except Exception:
            pass

    return start_time, end_time, duration_sec

def upsert_stats(out_file, row):
    out_file.parent.mkdir(parents=True, exist_ok=True)
    existing = []
    if out_file.exists():
        with out_file.open() as fh:
            existing = list(csv.DictReader(fh, delimiter="\t"))
    
    kept = existing[:]
    found = False
    for r in kept:
        if r.get("mode") == row["mode"]:
            r.update(row)
            found = True
            break
    
    if not found:
        kept.append(row)

    with out_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TRIM_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)

def main():
    valid_srrs = set(get_srrs())
    all_timings = []
    
    stats_written = 0
    
    # Use glob to find files quickly without deep traversing
    err_files = []
    for base_dir in [FB_RESULTS_DIR, FB_ORIGINAL_DIR]:
        if base_dir.exists():
            err_files.extend(list(base_dir.glob("*/*/_out_err/*.err")))
            err_files.extend(list(base_dir.glob("*/*/*trimmomatic/_out_err/*.err")))
            err_files.extend(list(base_dir.glob("*/*/*trimmomatic/*/*.err")))
            
    processed_stems = set()
    
    for err_file in err_files:
        if err_file.stem in processed_stems:
            continue
        
        # Extract project and srr from path
        # Example: .../PRJNA664293/SRR12676593/_out_err/run_trimming_...
        try:
            parts = err_file.parts
            out_err_idx = parts.index("_out_err")
            srr_id = parts[out_err_idx - 1]
            # Handle the case where there is an intermediate folder like <SRR>_trimmomatic
            if srr_id.endswith("_trimmomatic"):
                srr_id = parts[out_err_idx - 2]
                project_id = parts[out_err_idx - 3]
            else:
                project_id = parts[out_err_idx - 2]
        except ValueError:
            # Maybe not in _out_err
            continue
            
        if (project_id, srr_id) not in valid_srrs:
            continue
            
        processed_stems.add(err_file.stem)
        
        try:
            err_text = err_file.read_text(encoding="utf-8", errors="replace")
            stat_row = parse_trimmomatic_stats(err_text, project_id, srr_id)
            if not stat_row:
                continue
            
            # Write stat row
            stats_file = TRIM_BASE / project_id / f"{srr_id}_trimmomatic_stats.tsv"
            upsert_stats(stats_file, stat_row)
            stats_written += 1
            
            # Check for out file
            out_file = err_file.with_suffix(".out")
            if out_file.exists():
                out_text = out_file.read_text(encoding="utf-8", errors="replace")
                start_t, end_t, dur = parse_timing(out_text, stat_row["mode"])
                all_timings.append({
                    "project_id": project_id,
                    "srr_id": srr_id,
                    "mode": stat_row["mode"],
                    "start_time": start_t,
                    "end_time": end_t,
                    "duration_sec": dur
                })
                
        except Exception as e:
            print(f"Error processing {err_file}: {e}")

    # Write timings
    TIMING_FILE.parent.mkdir(parents=True, exist_ok=True)
    existing_timings = []
    if TIMING_FILE.exists():
        with TIMING_FILE.open() as fh:
            existing_timings = list(csv.DictReader(fh, delimiter="\t"))
    
    existing_keys = {(r.get("project_id", ""), r.get("srr_id", ""), r.get("mode", "")) for r in existing_timings}
    kept_timings = existing_timings[:]
    
    for row in all_timings:
        key = (row["project_id"], row["srr_id"], row["mode"])
        if key not in existing_keys:
            kept_timings.append(row)
            existing_keys.add(key)
    
    with TIMING_FILE.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIMING_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept_timings)
        
    print(f"Stats rows updated/written: {stats_written}")
    print(f"Timing rows total (including existing): {len(kept_timings)}")

if __name__ == "__main__":
    main()
