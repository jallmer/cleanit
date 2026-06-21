#!/usr/bin/env python3
"""
Merge ja and fb bowtie2 stats and extract fb timings.
"""

import sqlite3
import os
import re
import csv
from datetime import datetime
from pathlib import Path

DB_FILE = os.path.expanduser("~/srr_queue.db")
FB_RESULTS_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_trimming/records/fb_results")
FB_ORIGINAL_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/fb/omiks_project/results")
JA_STATS_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/flattened_bowtie2_stats")
FINAL_RECORDS = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_alignment/records")

BOWTIE_HEADER = [
    "project_id",
    "srr_id",
    "mode",
    "layout",
    "total_reads",
    "paired_reads",
    "unpaired_reads",
    "aligned_exactly_1",
    "aligned_gt1",
    "aligned_0",
    "concordant_exactly_1",
    "concordant_gt1",
    "discordant_1",
    "pairs_0_concordant",
    "overall_alignment_rate",
]

COMBINED_HEADER = ["source"] + BOWTIE_HEADER
TIMING_HEADER = ["project_id", "srr_id", "mode", "start_time", "end_time", "duration_sec"]

def get_srrs():
    conn = sqlite3.connect(DB_FILE)
    cur = conn.cursor()
    cur.execute("SELECT project_id, srr_id FROM srr_queue")
    srrs = cur.fetchall()
    conn.close()
    return srrs

def parse_timing(out_text):
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

def parse_bowtie_err(err_text, project_id, srr_id, mode):
    total_reads = re.search(r"(\d+)\s+reads; of these:", err_text, re.M)
    overall = re.search(r"([\d.]+%)\s+overall alignment rate", err_text, re.M)
    if not total_reads or not overall:
        return None

    layout = "SE" if "were unpaired;" in err_text else "PE"

    row = {
        "project_id": project_id,
        "srr_id": srr_id,
        "mode": mode,
        "layout": layout,
        "total_reads": total_reads.group(1),
        "paired_reads": "",
        "unpaired_reads": "",
        "aligned_exactly_1": "",
        "aligned_gt1": "",
        "aligned_0": "",
        "concordant_exactly_1": "",
        "concordant_gt1": "",
        "discordant_1": "",
        "pairs_0_concordant": "",
        "overall_alignment_rate": overall.group(1),
    }

    if layout == "PE":
        patterns = {
            "paired_reads": r"(\d+)\s+\([^)]+\)\s+were paired; of these:",
            "unpaired_reads": r"(\d+)\s+\([^)]+\)\s+were unpaired; of these:",
            "concordant_exactly_1": r"(\d+)\s+\([^)]+\)\s+aligned concordantly exactly 1 time",
            "concordant_gt1": r"(\d+)\s+\([^)]+\)\s+aligned concordantly >1 times",
            "discordant_1": r"(\d+)\s+\([^)]+\)\s+aligned discordantly 1 time",
            "pairs_0_concordant": r"(\d+)(?:\s+\([^)]+\))?\s+pairs aligned concordantly 0 times",
        }
    else:
        patterns = {
            "aligned_exactly_1": r"(\d+)\s+\([^)]+\)\s+aligned exactly 1 time$",
            "aligned_gt1": r"(\d+)\s+\([^)]+\)\s+aligned >1 times$",
            "aligned_0": r"(\d+)\s+\([^)]+\)\s+aligned 0 times$",
        }

    for key, pattern in patterns.items():
        m = re.search(pattern, err_text, re.M)
        if m:
            row[key] = m.group(1)

    return row

def consolidate_ja_stats():
    ja_rows = []
    for f in JA_STATS_DIR.glob("*/*_bowtie2_stats.tsv"):
        with f.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                # Add source
                row["source"] = "ja"
                ja_rows.append(row)
    return ja_rows

def main():
    FINAL_RECORDS.mkdir(parents=True, exist_ok=True)
    valid_srrs = set(get_srrs())
    
    ja_rows = consolidate_ja_stats()
    
    fb_rows = []
    fb_timings = []
    
    # Process fb logs
    err_files = []
    for base_dir in [FB_RESULTS_DIR, FB_ORIGINAL_DIR]:
        if base_dir.exists():
            err_files.extend(list(base_dir.glob("*/*/_out_err/run_*alignment*.err")))
            err_files.extend(list(base_dir.glob("*/*/*trimmomatic/_out_err/run_*alignment*.err")))
            err_files.extend(list(base_dir.glob("*/*/*trimmomatic/*/*.err")))

    processed_stems = set()
    
    for err_file in err_files:
        if err_file.stem in processed_stems:
            continue
            
        name = err_file.name
        if not ("alignment" in name):
            continue

        try:
            parts = err_file.parts
            out_err_idx = parts.index("_out_err")
            srr_id = parts[out_err_idx - 1]
            if srr_id.endswith("_trimmomatic"):
                srr_id = parts[out_err_idx - 2]
                project_id = parts[out_err_idx - 3]
            else:
                project_id = parts[out_err_idx - 2]
        except ValueError:
            continue

        if (project_id, srr_id) not in valid_srrs:
            continue

        processed_stems.add(err_file.stem)
        
        err_text = err_file.read_text(encoding="utf-8", errors="replace")
        out_file = err_file.with_suffix(".out")
        out_text = ""
        if out_file.exists():
            out_text = out_file.read_text(encoding="utf-8", errors="replace")
        
        mode = "unknown"
        if out_text:
            if re.search(r"aligned .*/raw_data/", out_text):
                mode = "untrimmed"
            else:
                m = re.search(r"aligned .*trimmomatic_([a-zA-Z0-9]+)\s+to", out_text)
                if m:
                    mode = m.group(1)
                    if mode == "adapter":
                        mode = "adapter_only"

        # Fallback if the out file is missing or doesn't have the aligned line
        if mode == "unknown":
            if "trmd_alignment" in name:
                if "P35" in str(err_file): mode = "P35"
                elif "P20" in str(err_file): mode = "P20"
                elif "P10" in str(err_file): mode = "P10"
                elif "P5" in str(err_file): mode = "P5"
                elif "adapter" in str(err_file): mode = "adapter_only"
            else:
                mode = "untrimmed"

        stat_row = parse_bowtie_err(err_text, project_id, srr_id, mode)
        if stat_row:
            fb_stats_row = stat_row.copy()
            fb_stats_row["source"] = "fb"
            fb_rows.append(fb_stats_row)
            
            # Timing
            if out_text:
                start_t, end_t, dur = parse_timing(out_text)
                if dur:
                    fb_timings.append({
                        "project_id": project_id,
                        "srr_id": srr_id,
                        "mode": mode,
                        "start_time": start_t,
                        "end_time": end_t,
                        "duration_sec": dur
                    })

    # Write JA table
    with (FINAL_RECORDS / "ja_alignment_stats.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["source"] + BOWTIE_HEADER, delimiter="\t", extrasaction='ignore')
        writer.writeheader()
        writer.writerows(ja_rows)

    # Write FB table
    with (FINAL_RECORDS / "fb_alignment_stats.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["source"] + BOWTIE_HEADER, delimiter="\t", extrasaction='ignore')
        writer.writeheader()
        writer.writerows(fb_rows)

    # Write Combined table
    with (FINAL_RECORDS / "combined_alignment_stats.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COMBINED_HEADER, delimiter="\t", extrasaction='ignore')
        writer.writeheader()
        writer.writerows(ja_rows + fb_rows)

    # Write FB Timing table
    with (FINAL_RECORDS / "fb_alignment_timing.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIMING_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(fb_timings)

    print(f"JA stats rows: {len(ja_rows)}")
    print(f"FB stats rows: {len(fb_rows)}")
    print(f"FB timing rows: {len(fb_timings)}")
    print(f"Combined stats rows: {len(ja_rows) + len(fb_rows)}")

if __name__ == "__main__":
    main()
