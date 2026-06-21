#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

root = Path("/scratch/hpc-prf-omiks/ja/aligner_benchmark")
speed = root / "results" / "aligner_speed.tsv"
trim = root / "results" / "trimming_speed.tsv"
out = root / "results" / "aligner_speed_summary.tsv"

if not speed.exists():
    raise SystemExit(f"Missing {speed}")

df = pd.read_csv(speed, sep="\t")
ok = df[df["status"].eq("ok")].copy()
if ok.empty:
    raise SystemExit("No completed aligner rows yet")

summary = (
    ok.groupby(["aligner", "mode"], dropna=False)
    .agg(
        n=("srr_id", "nunique"),
        median_duration_sec=("duration_sec", "median"),
        mean_duration_sec=("duration_sec", "mean"),
        min_duration_sec=("duration_sec", "min"),
        max_duration_sec=("duration_sec", "max"),
        median_reads_M=("total_reads_M", "median"),
    )
    .reset_index()
)
summary.to_csv(out, sep="\t", index=False)

if trim.exists():
    trim_df = pd.read_csv(trim, sep="\t")
    trim_df.to_csv(root / "results" / "trimming_speed_copy.tsv", sep="\t", index=False)

print(f"Wrote {out}")
