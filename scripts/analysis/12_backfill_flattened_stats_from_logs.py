#!/usr/bin/env python3
"""
Backfill flattened per-SRR Bowtie2, Trimmomatic, and timing stats from existing
Slurm logs.

Reads:
  ~/logs/v2_pipeline_*.out
  ~/logs/v2_pipeline_*.err

Writes:
  /scratch/hpc-prf-omiks/ja/flattened_trimmomatic_stats/<PROJECT>/<SRR>_trimmomatic_stats.tsv
  /scratch/hpc-prf-omiks/ja/flattened_bowtie2_stats/<PROJECT>/<SRR>_bowtie2_stats.tsv
  /scratch/hpc-prf-omiks/ja/flattened_timings/<PROJECT>/<SRR>_timings.tsv

Behavior:
  - skips active running v2_pipeline job IDs
  - writes one row per trimming mode
  - preserves existing flattened rows; only fills missing SRR/mode or stage rows
  - optionally deletes completed logs after successful extraction
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import subprocess
from pathlib import Path


HOME = Path.home()
LOG_DIR = HOME / "logs"
TRIM_BASE = Path("/scratch/hpc-prf-omiks/ja/flattened_trimmomatic_stats")
BOWTIE_BASE = Path("/scratch/hpc-prf-omiks/ja/flattened_bowtie2_stats")
TIMING_BASE = Path("/scratch/hpc-prf-omiks/ja/flattened_timings")
MODE_ORDER = ["adapter_only", "P5", "P10", "P20", "P35"]

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

TIMING_HEADER = ["stage", "start_epoch", "end_epoch", "duration_sec"]


def running_pipeline_job_ids() -> set[str]:
    try:
        result = subprocess.run(
            ["bash", "-lc", "squeue -u $(whoami) -h -t R -n v2_pipeline -o '%i'"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return set()
    return {line.strip() for line in result.stdout.splitlines() if line.strip()}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def parse_out_metadata(out_text: str) -> tuple[str | None, str | None, str]:
    m = re.search(r"Target SRR:\s+(\S+)\s+\(([^)]+)\)", out_text)
    if not m:
        return None, None, "SE"
    srr_id, project_id = m.group(1), m.group(2)
    layout = "PE" if "Parsed Layout: PE" in out_text else "SE"
    return project_id, srr_id, layout


def extract_trim_rows(err_text: str, project_id: str, srr_id: str, layout: str) -> list[dict[str, str]]:
    rows = []
    pe_re = re.compile(
        r"Input Read Pairs:\s+(\d+)\s+Both Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+"
        r"Forward Only Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+"
        r"Reverse Only Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Dropped:\s+(\d+)\s+\(([\d.]+)%\)"
    )
    se_re = re.compile(
        r"Input Reads:\s+(\d+)\s+Surviving:\s+(\d+)\s+\(([\d.]+)%\)\s+Dropped:\s+(\d+)\s+\(([\d.]+)%\)"
    )

    if layout == "PE":
        matches = pe_re.findall(err_text)
        for mode, groups in zip(MODE_ORDER, matches):
            input_reads = int(groups[0])
            both = int(groups[1])
            forward = int(groups[3])
            reverse = int(groups[5])
            dropped = int(groups[7])
            surviving = both + forward + reverse
            surviving_pct = f"{(surviving / input_reads * 100) if input_reads else 0:.2f}"
            rows.append(
                {
                    "project_id": project_id,
                    "srr_id": srr_id,
                    "mode": mode,
                    "layout": layout,
                    "input_reads": str(input_reads),
                    "surviving": str(surviving),
                    "surviving_pct": surviving_pct,
                    "dropped": str(dropped),
                    "dropped_pct": groups[8],
                }
            )
    else:
        matches = se_re.findall(err_text)
        for mode, groups in zip(MODE_ORDER, matches):
            rows.append(
                {
                    "project_id": project_id,
                    "srr_id": srr_id,
                    "mode": mode,
                    "layout": layout,
                    "input_reads": groups[0],
                    "surviving": groups[1],
                    "surviving_pct": groups[2],
                    "dropped": groups[3],
                    "dropped_pct": groups[4],
                }
            )

    return rows


def extract_bowtie_rows(err_text: str, project_id: str, srr_id: str, layout: str) -> list[dict[str, str]]:
    rows = []
    blocks = re.split(r"Executing Bowtie2 Matrix Arrays Mapping to Samtools Index\.\.\.\n", err_text)
    for mode, block in zip(["untrimmed"] + MODE_ORDER, blocks[1:]):
        total_reads = re.search(r"(\d+)\s+reads; of these:", block, re.M)
        overall = re.search(r"([\d.]+%)\s+overall alignment rate", block, re.M)
        if not total_reads or not overall:
            continue

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
            m = re.search(pattern, block, re.M)
            if m:
                row[key] = m.group(1)

        rows.append(row)

    return rows


def extract_timing_rows(out_text: str) -> list[dict[str, str]]:
    rows = []
    seen = set()
    for name, duration in re.findall(r"^\[TIMING\]\s+([A-Za-z0-9_]+):\s+(\d+)s\s*$", out_text, re.M):
        # Log-derived timings do not preserve the relative start/end counters.
        # Downstream throughput code only requires duration_sec.
        if name in seen:
            continue
        seen.add(name)
        rows.append(
            {
                "stage": name,
                "start_epoch": "",
                "end_epoch": "",
                "duration_sec": duration,
            }
        )

    total = re.search(r"PIPELINE COMPLETED SUCCESSFULLY \(total:\s*(\d+)s\)", out_text)
    if total and "pipeline_total" not in seen:
        rows.append(
            {
                "stage": "pipeline_total",
                "start_epoch": "",
                "end_epoch": "",
                "duration_sec": total.group(1),
            }
        )
    return rows


def upsert_rows(out_file: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    existing: list[dict[str, str]] = []
    if out_file.exists():
        with out_file.open() as fh:
            existing = list(csv.DictReader(fh, delimiter="\t"))

    existing_keys = {
        (row.get("project_id", ""), row.get("srr_id", ""), row.get("mode", ""))
        for row in existing
    }
    kept = existing[:]
    for row in existing:
        for key in header:
            row.setdefault(key, "")
    for row in rows:
        key = (row["project_id"], row["srr_id"], row["mode"])
        if key not in existing_keys:
            kept.append(row)
            existing_keys.add(key)

    with out_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)


def upsert_timing_rows(out_file: Path, rows: list[dict[str, str]]) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    existing: list[dict[str, str]] = []
    if out_file.exists():
        with out_file.open() as fh:
            existing = list(csv.DictReader(fh, delimiter="\t"))

    existing_stages = {row.get("stage", "") for row in existing}
    kept = existing[:]
    for row in existing:
        for key in TIMING_HEADER:
            row.setdefault(key, "")
    for row in rows:
        stage = row["stage"]
        if stage not in existing_stages:
            kept.append(row)
            existing_stages.add(stage)

    with out_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIMING_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)


def is_completed(out_text: str) -> bool:
    return "PIPELINE COMPLETED SUCCESSFULLY" in out_text


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delete-extracted", action="store_true")
    args = parser.parse_args()

    active = running_pipeline_job_ids()
    err_files = sorted(Path(p) for p in glob.glob(str(LOG_DIR / "v2_pipeline_*.err")))

    processed = 0
    deleted = 0
    skipped_active = 0
    skipped_missing_meta = 0
    trim_files_written = 0
    bowtie_files_written = 0
    timing_files_written = 0

    for err_file in err_files:
        job_id = err_file.stem.replace("v2_pipeline_", "")
        if job_id in active:
            skipped_active += 1
            continue

        out_file = err_file.with_suffix(".out")
        if not out_file.exists():
            continue

        out_text = read_text(out_file)
        err_text = read_text(err_file)
        project_id, srr_id, layout = parse_out_metadata(out_text)
        if not project_id or not srr_id:
            skipped_missing_meta += 1
            continue

        trim_rows = extract_trim_rows(err_text, project_id, srr_id, layout)
        bowtie_rows = extract_bowtie_rows(err_text, project_id, srr_id, layout)
        timing_rows = extract_timing_rows(out_text)

        if trim_rows:
            upsert_rows(
                TRIM_BASE / project_id / f"{srr_id}_trimmomatic_stats.tsv",
                TRIM_HEADER,
                trim_rows,
            )
            trim_files_written += 1
        if bowtie_rows:
            upsert_rows(
                BOWTIE_BASE / project_id / f"{srr_id}_bowtie2_stats.tsv",
                BOWTIE_HEADER,
                bowtie_rows,
            )
            bowtie_files_written += 1
        if timing_rows:
            upsert_timing_rows(
                TIMING_BASE / project_id / f"{srr_id}_timings.tsv",
                timing_rows,
            )
            timing_files_written += 1

        processed += 1

        if args.delete_extracted and is_completed(out_text):
            err_file.unlink(missing_ok=True)
            out_file.unlink(missing_ok=True)
            deleted += 1

    print(f"processed_logs\t{processed}")
    print(f"trim_files_written\t{trim_files_written}")
    print(f"bowtie_files_written\t{bowtie_files_written}")
    print(f"timing_files_written\t{timing_files_written}")
    print(f"skipped_active\t{skipped_active}")
    print(f"skipped_missing_meta\t{skipped_missing_meta}")
    print(f"deleted_completed_logs\t{deleted}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
