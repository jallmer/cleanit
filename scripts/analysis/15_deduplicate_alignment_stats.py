#!/usr/bin/env python3
"""
Deduplicate Bowtie2 alignment stats and timings.
"""

import pandas as pd
from pathlib import Path

FINAL_RECORDS = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_alignment/records")

def main():
    # 1. Load fb datasets
    fb_stats = pd.read_csv(FINAL_RECORDS / "fb_alignment_stats.tsv", sep="\t")
    fb_timings = pd.read_csv(FINAL_RECORDS / "fb_alignment_timing.tsv", sep="\t")
    
    assert len(fb_stats) == len(fb_timings), "Length mismatch between fb stats and timings"
    
    # 2. Add timing columns to stats for deduplication
    fb_stats["start_time"] = fb_timings["start_time"]
    fb_stats["end_time_str"] = fb_timings["end_time"]
    fb_stats["duration_sec"] = fb_timings["duration_sec"]
    
    # Parse for sorting
    fb_stats["end_time_dt"] = pd.to_datetime(fb_timings["end_time"], format="%Y-%m-%dT%H:%M:%S%z", errors="coerce")
    
    # 3. Sort by end_time descending (so latest is first)
    fb_stats = fb_stats.sort_values(by="end_time_dt", ascending=False)
    
    # 4. Drop duplicates by project, srr, mode, keeping the first (which is the latest due to sort)
    fb_stats_dedup = fb_stats.drop_duplicates(subset=["project_id", "srr_id", "mode"], keep="first")
    
    # 5. Extract and rename timing columns
    fb_timings_dedup = fb_stats_dedup[["project_id", "srr_id", "mode", "start_time", "end_time_str", "duration_sec"]].copy()
    fb_timings_dedup.rename(columns={"end_time_str": "end_time"}, inplace=True)
    
    # 6. Clean stats dataframe
    stats_cols = [c for c in fb_stats_dedup.columns if c not in ["start_time", "end_time_str", "duration_sec", "end_time_dt"]]
    fb_stats_dedup = fb_stats_dedup[stats_cols]
    
    # 7. Load ja stats and deduplicate within it
    ja_stats = pd.read_csv(FINAL_RECORDS / "ja_alignment_stats.tsv", sep="\t")
    ja_stats = ja_stats.drop_duplicates(subset=["project_id", "srr_id", "mode"], keep="first")
    
    # 8. Merge ja and fb, giving precedence to ja
    # Concat with ja FIRST, then drop duplicates keeping "first"
    combined = pd.concat([ja_stats, fb_stats_dedup], ignore_index=True)
    combined_dedup = combined.drop_duplicates(subset=["project_id", "srr_id", "mode"], keep="first")
    
    # 9. Save out the clean files
    fb_stats_dedup.to_csv(FINAL_RECORDS / "fb_alignment_stats.tsv", sep="\t", index=False)
    fb_timings_dedup.to_csv(FINAL_RECORDS / "fb_alignment_timing.tsv", sep="\t", index=False)
    ja_stats.to_csv(FINAL_RECORDS / "ja_alignment_stats.tsv", sep="\t", index=False)
    combined_dedup.to_csv(FINAL_RECORDS / "combined_alignment_stats.tsv", sep="\t", index=False)
    
    print("Deduplication complete:")
    print(f"  JA deduplicated rows: {len(ja_stats)}")
    print(f"  FB deduplicated rows: {len(fb_stats_dedup)}")
    print(f"  FB deduplicated timing rows: {len(fb_timings_dedup)}")
    print(f"  Combined deduplicated rows: {len(combined_dedup)}")

if __name__ == '__main__':
    main()
