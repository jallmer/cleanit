#!/usr/bin/env python3
"""Build a consolidated FastQC bundle from fb and ja flattened outputs.

The bundle intentionally excludes sequence/alignment payloads. It keeps FastQC
reports, extracted FastQC report files, FastQC/get-data logs, and timing
summaries that can be parsed from existing logs/timing tables.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


BASE = Path("/pc2/users/o/omiks001/hpc-prf-omiks")
FB_RESULTS = BASE / "fb/omiks_project/results"
JA = BASE / "ja"
FLATTENED_FASTQC = JA / "flattened_fastqc_raw"
FLATTENED_TIMINGS = JA / "flattened_timings"
OUT = JA / "final_fastqc"

FORBIDDEN_EXTENSIONS = (
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".sam",
    ".sam.gz",
    ".bam",
    ".bam.gz",
    ".cram",
    ".cram.gz",
    ".bai",
    ".crai",
)


def is_forbidden_payload(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in FORBIDDEN_EXTENSIONS)


def safe_copy(src: Path, dst: Path) -> tuple[str, int]:
    """Copy src to dst, returning action and previous file size if overwritten."""
    if is_forbidden_payload(src):
        raise ValueError(f"Refusing to copy payload file: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    overwritten_size = dst.stat().st_size if dst.exists() else -1
    shutil.copy2(src, dst)
    return ("overwrote" if overwritten_size >= 0 else "created", overwritten_size)


def fb_fastqc_files() -> list[Path]:
    if not FB_RESULTS.exists():
        return []
    cmd = [
        "find",
        str(FB_RESULTS),
        "-type",
        "f",
        "(",
        "-iname",
        "*_fastqc.zip",
        "-o",
        "-iname",
        "*_fastqc.html",
        "-o",
        "-path",
        "*_fastqc/*",
        "-o",
        "-path",
        "*/_out_err/run_get_SRR_data_and_FastQc*",
        ")",
    ]
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
    return sorted(Path(line) for line in proc.stdout.splitlines() if line and not is_forbidden_payload(Path(line)))


def flattened_fastqc_files() -> list[Path]:
    if not FLATTENED_FASTQC.exists():
        return []
    cmd = ["find", str(FLATTENED_FASTQC), "-type", "f"]
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
    return sorted(Path(line) for line in proc.stdout.splitlines() if line and not is_forbidden_payload(Path(line)))


def project_srr_from_fb(path: Path) -> tuple[str, str]:
    rel = path.relative_to(FB_RESULTS)
    project = rel.parts[0] if len(rel.parts) > 0 else ""
    srr = rel.parts[1] if len(rel.parts) > 1 else ""
    return project, srr


def project_srr_from_flattened(path: Path) -> tuple[str, str]:
    rel = path.relative_to(FLATTENED_FASTQC)
    project = rel.parts[0] if len(rel.parts) > 0 else ""
    match = re.search(r"(SRR\d+)", "/".join(rel.parts))
    return project, match.group(1) if match else ""


def parse_isoish(text: str) -> datetime | None:
    for fmt in ("%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%dT%H:%M:%S"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            pass
    return None


def parse_fb_logs() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if not FB_RESULTS.exists():
        return rows
    pattern = re.compile(r"run_get_SRR_data_and_FastQc_\d+\.out$")
    for path in sorted(FB_RESULTS.rglob("run_get_SRR_data_and_FastQc_*.out")):
        if not pattern.search(path.name):
            continue
        project, srr = project_srr_from_fb(path)
        text = path.read_text(errors="replace")
        start = None
        end = None
        download_start = None
        download_end = None
        skipped_download = False
        skipped_dumping = False
        ran_fastqc = False
        ended_cleanly = False
        err_path = path.with_suffix(".err")
        err_text = err_path.read_text(errors="replace") if err_path.exists() else ""
        fastqc_complete_count = text.count("Analysis complete for ") + err_text.count("Analysis complete for ")
        for line in text.splitlines():
            if "Startet get_SRR_data_and_FastQc.sh" in line:
                m = re.search(r" at (\S+)$", line)
                if m:
                    start = parse_isoish(m.group(1))
            elif "ended get_SRR_data_and_FastQc.sh" in line:
                m = re.search(r" at (\S+)$", line)
                if m:
                    end = parse_isoish(m.group(1))
            elif "Downloading via HTTPS" in line and download_start is None:
                m = re.match(r"(\S+)", line)
                if m:
                    download_start = parse_isoish(m.group(1))
            elif "HTTPS download succeed" in line:
                m = re.match(r"(\S+)", line)
                if m:
                    download_end = parse_isoish(m.group(1))
            elif "already fully exist" in line and "Skipping download" in line:
                skipped_download = True
            elif "Skipping dumping" in line:
                skipped_dumping = True
            elif "Performing FastQc" in line:
                ran_fastqc = True
            elif "ended get_SRR_data_and_FastQc.sh" in line:
                ended_cleanly = True
        job_wall_sec = int((end - start).total_seconds()) if start and end else ""
        fastqc_only_approx_sec = (
            job_wall_sec
            if job_wall_sec != ""
            and skipped_download
            and skipped_dumping
            and ran_fastqc
            and ended_cleanly
            and fastqc_complete_count > 0
            else ""
        )
        rows.append(
            {
                "project_id": project,
                "srr_id": srr,
                "log_path": str(path),
                "job_wall_sec": job_wall_sec,
                "sra_https_download_sec": int((download_end - download_start).total_seconds())
                if download_start and download_end
                else "",
                "skipped_download": "yes" if skipped_download else "no",
                "skipped_dumping": "yes" if skipped_dumping else "no",
                "ran_fastqc": "yes" if ran_fastqc else "no",
                "fastqc_complete_count": fastqc_complete_count,
                "fb_cpus_per_task": 1,
                "fastqc_only_approx_sec": fastqc_only_approx_sec,
                "fastqc_only_approx_basis": "whole fb job wall time; download and dumping skipped; FastQC completed"
                if fastqc_only_approx_sec != ""
                else "",
            }
        )
    return rows


def parse_flattened_fastqc_timings() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if not FLATTENED_TIMINGS.exists():
        return rows
    for path in sorted(FLATTENED_TIMINGS.glob("*/*_timings.tsv")):
        project = path.parent.name
        srr = path.name.removesuffix("_timings.tsv")
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row.get("stage") == "fastqc":
                    rows.append(
                        {
                            "project_id": project,
                            "srr_id": srr,
                            "timing_path": str(path),
                            "start_epoch": row.get("start_epoch", ""),
                            "end_epoch": row.get("end_epoch", ""),
                            "duration_sec": row.get("duration_sec", ""),
                        }
                    )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []

    fb_files = fb_fastqc_files()
    for src in fb_files:
        rel = src.relative_to(FB_RESULTS)
        dst = OUT / "reports" / rel
        action, previous_size = safe_copy(src, dst)
        project, srr = project_srr_from_fb(src)
        manifest.append(
            {
                "source_priority": 1,
                "source": "fb_results",
                "source_kind": "fastqc_or_get_data_log",
                "project_id": project,
                "srr_id": srr,
                "source_path": str(src),
                "destination_path": str(dst),
                "action": action,
                "source_size_bytes": src.stat().st_size,
                "previous_destination_size_bytes": previous_size if previous_size >= 0 else "",
            }
        )

    flat_files = flattened_fastqc_files()
    for src in flat_files:
        rel = src.relative_to(FLATTENED_FASTQC)
        dst = OUT / "reports" / rel
        action, previous_size = safe_copy(src, dst)
        project, srr = project_srr_from_flattened(src)
        manifest.append(
            {
                "source_priority": 2,
                "source": "flattened_fastqc_raw",
                "source_kind": "fastqc_report",
                "project_id": project,
                "srr_id": srr,
                "source_path": str(src),
                "destination_path": str(dst),
                "action": action,
                "source_size_bytes": src.stat().st_size,
                "previous_destination_size_bytes": previous_size if previous_size >= 0 else "",
            }
        )

    manifest_fields = [
        "source_priority",
        "source",
        "source_kind",
        "project_id",
        "srr_id",
        "source_path",
        "destination_path",
        "action",
        "source_size_bytes",
        "previous_destination_size_bytes",
    ]
    write_tsv(OUT / "manifest.tsv", manifest, manifest_fields)

    log_rows = parse_fb_logs()
    write_tsv(
        OUT / "timing/fb_get_srr_fastqc_log_summary.tsv",
        log_rows,
        [
            "project_id",
            "srr_id",
            "log_path",
            "job_wall_sec",
            "sra_https_download_sec",
            "skipped_download",
            "skipped_dumping",
            "ran_fastqc",
            "fastqc_complete_count",
            "fb_cpus_per_task",
            "fastqc_only_approx_sec",
            "fastqc_only_approx_basis",
        ],
    )

    timing_rows = parse_flattened_fastqc_timings()
    write_tsv(
        OUT / "timing/flattened_fastqc_stage_timings.tsv",
        timing_rows,
        ["project_id", "srr_id", "timing_path", "start_epoch", "end_epoch", "duration_sec"],
    )

    summary = [
        {"metric": "fb_fastqc_related_files_found", "value": len(fb_files)},
        {"metric": "flattened_fastqc_files_found", "value": len(flat_files)},
        {"metric": "manifest_rows", "value": len(manifest)},
        {"metric": "manifest_overwrites_by_flattened_layer", "value": sum(1 for r in manifest if r["source"] == "flattened_fastqc_raw" and r["action"] == "overwrote")},
        {"metric": "fb_get_srr_fastqc_logs_with_wall_time", "value": sum(1 for r in log_rows if r["job_wall_sec"] != "")},
        {"metric": "fb_get_srr_fastqc_logs_with_sra_download_time", "value": sum(1 for r in log_rows if r["sra_https_download_sec"] != "")},
        {"metric": "flattened_fastqc_stage_timing_rows", "value": len(timing_rows)},
    ]
    write_tsv(OUT / "summary.tsv", summary, ["metric", "value"])

    print(f"final_fastqc: {OUT}")
    for row in summary:
        print(f"{row['metric']}\t{row['value']}")


if __name__ == "__main__":
    main()
