#!/usr/bin/env python3
"""Extract one row of crucial FastQC/timing metadata per DB SRR with a report."""

from __future__ import annotations

import csv
import math
import re
import sqlite3
from collections import defaultdict
from pathlib import Path


ROOT = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja")
FINAL = ROOT / "final_fastqc"
MANIFEST = FINAL / "manifest.tsv"
DB = Path("/pc2/users/o/omiks001/srr_queue.db")
OUT = FINAL / "srr_fastqc_metrics.tsv"


def srr_from_text(text: str) -> str:
    m = re.search(r"(SRR\d+)", text)
    return m.group(1) if m else ""


def read_mate_from_name(text: str) -> str:
    if re.search(r"_1_fastqc", text):
        return "R1"
    if re.search(r"_2_fastqc", text):
        return "R2"
    return "single_or_unlabelled"


def as_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def median(values: list[float]) -> float | None:
    values = sorted(v for v in values if v is not None and not math.isnan(v))
    if not values:
        return None
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2


def load_db() -> dict[str, dict[str, str]]:
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(
            """
            select project_id, srr_id, status, run_status, fetch_lane, size_gb,
                   sra_size_bytes, fastq_size_bytes
            from srr_queue
            """
        ).fetchall()
        return {r["srr_id"]: dict(r) for r in rows}
    finally:
        con.close()


def parse_fastqc(path: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "fastqc_version": "",
        "filename": "",
        "file_type": "",
        "encoding": "",
        "total_sequences": "",
        "total_bases": "",
        "sequences_flagged_poor_quality": "",
        "sequence_length": "",
        "gc_percent": "",
        "per_base_quality_status": "",
        "adapter_content_status": "",
        "sequence_duplication_status": "",
        "overrepresented_sequences_status": "",
        "mean_quality_all_bases": "",
        "mean_quality_first_bin": "",
        "mean_quality_terminal_bin": "",
        "mean_quality_last_10_bins": "",
        "median_quality_terminal_bin": "",
        "lowest_mean_quality_bin": "",
        "lowest_mean_quality_position": "",
    }
    lines = path.read_text(errors="replace").splitlines()
    if lines and lines[0].startswith("##FastQC"):
        parts = lines[0].split("\t")
        result["fastqc_version"] = parts[1] if len(parts) > 1 else ""

    module = ""
    basic: dict[str, str] = {}
    quality_rows: list[tuple[str, float, float | None]] = []

    for line in lines:
        if line.startswith(">>") and not line.startswith(">>END_MODULE"):
            bits = line[2:].split("\t")
            module = bits[0]
            status = bits[1] if len(bits) > 1 else ""
            if module == "Per base sequence quality":
                result["per_base_quality_status"] = status
            elif module == "Adapter Content":
                result["adapter_content_status"] = status
            elif module == "Sequence Duplication Levels":
                result["sequence_duplication_status"] = status
            elif module == "Overrepresented sequences":
                result["overrepresented_sequences_status"] = status
            continue
        if line.startswith(">>END_MODULE"):
            module = ""
            continue
        if not line or line.startswith("#"):
            continue
        if module == "Basic Statistics":
            parts = line.split("\t")
            if len(parts) >= 2:
                basic[parts[0]] = parts[1]
        elif module == "Per base sequence quality":
            parts = line.split("\t")
            if len(parts) >= 3:
                mean = as_float(parts[1])
                med = as_float(parts[2])
                if mean is not None:
                    quality_rows.append((parts[0], mean, med))

    result["filename"] = basic.get("Filename", "")
    result["file_type"] = basic.get("File type", "")
    result["encoding"] = basic.get("Encoding", "")
    result["total_sequences"] = basic.get("Total Sequences", "")
    result["total_bases"] = basic.get("Total Bases", "")
    result["sequences_flagged_poor_quality"] = basic.get("Sequences flagged as poor quality", "")
    result["sequence_length"] = basic.get("Sequence length", "")
    result["gc_percent"] = basic.get("%GC", "")

    if quality_rows:
        means = [r[1] for r in quality_rows]
        result["mean_quality_all_bases"] = round(sum(means) / len(means), 4)
        result["mean_quality_first_bin"] = round(quality_rows[0][1], 4)
        result["mean_quality_terminal_bin"] = round(quality_rows[-1][1], 4)
        result["mean_quality_last_10_bins"] = round(sum(means[-10:]) / len(means[-10:]), 4)
        result["median_quality_terminal_bin"] = quality_rows[-1][2] if quality_rows[-1][2] is not None else ""
        low = min(quality_rows, key=lambda x: x[1])
        result["lowest_mean_quality_bin"] = round(low[1], 4)
        result["lowest_mean_quality_position"] = low[0]

    return result


def load_manifest_report_rows() -> dict[str, dict[str, str]]:
    """Return one retained fastqc_data row per SRR."""
    by_srr: dict[str, dict[str, str]] = {}
    with MANIFEST.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dst = row["destination_path"]
            if not dst.endswith("/fastqc_data.txt"):
                continue
            srr = row.get("srr_id") or srr_from_text(dst)
            by_srr[srr] = row
    return by_srr


def load_flattened_fastqc_timings() -> dict[str, dict[str, object]]:
    timings: dict[str, list[float]] = defaultdict(list)
    paths: dict[str, list[str]] = defaultdict(list)
    path = FINAL / "timing/flattened_fastqc_stage_timings.tsv"
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            val = as_float(row.get("duration_sec", ""))
            if val is not None:
                timings[row["srr_id"]].append(val)
                paths[row["srr_id"]].append(row.get("timing_path", ""))
    return {
        srr: {
            "fastqc_wall_sec_flattened_min": min(vals),
            "fastqc_wall_sec_flattened_max": max(vals),
            "fastqc_wall_sec_flattened_median": median(vals),
            "fastqc_timing_rows": len(vals),
            "fastqc_timing_path_example": paths[srr][0] if paths[srr] else "",
        }
        for srr, vals in timings.items()
    }


def load_fb_timing_summary() -> dict[str, dict[str, object]]:
    wall: dict[str, list[float]] = defaultdict(list)
    dl: dict[str, list[float]] = defaultdict(list)
    fastqc_approx: dict[str, list[float]] = defaultdict(list)
    log_paths: dict[str, list[str]] = defaultdict(list)
    path = FINAL / "timing/fb_get_srr_fastqc_log_summary.tsv"
    if not path.exists():
        return {}
    enhanced_rows: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            srr = row["srr_id"]
            w = as_float(row.get("job_wall_sec", ""))
            d = as_float(row.get("sra_https_download_sec", ""))
            if w is not None:
                wall[srr].append(w)
            if d is not None:
                dl[srr].append(d)
            log_paths[srr].append(row.get("log_path", ""))
            log_path = Path(row.get("log_path", ""))
            skipped_download = ""
            skipped_dumping = ""
            ran_fastqc = ""
            fastqc_complete_count = ""
            approx = None
            if log_path.exists():
                text = log_path.read_text(errors="replace")
                err_path = log_path.with_suffix(".err")
                err_text = err_path.read_text(errors="replace") if err_path.exists() else ""
                skipped_download_bool = "already fully exist" in text and "Skipping download" in text
                skipped_dumping_bool = "Skipping dumping" in text
                ran_fastqc_bool = "Performing FastQc" in text
                fastqc_complete_count_int = text.count("Analysis complete for ") + err_text.count("Analysis complete for ")
                ended_cleanly_bool = "ended get_SRR_data_and_FastQc.sh" in text
                skipped_download = "yes" if skipped_download_bool else "no"
                skipped_dumping = "yes" if skipped_dumping_bool else "no"
                ran_fastqc = "yes" if ran_fastqc_bool else "no"
                fastqc_complete_count = fastqc_complete_count_int
                if (
                    w is not None
                    and skipped_download_bool
                    and skipped_dumping_bool
                    and ran_fastqc_bool
                    and ended_cleanly_bool
                    and fastqc_complete_count_int > 0
                ):
                    approx = w
                    fastqc_approx[srr].append(w)
            enhanced = dict(row)
            enhanced.update(
                {
                    "skipped_download": skipped_download,
                    "skipped_dumping": skipped_dumping,
                    "ran_fastqc": ran_fastqc,
                    "fastqc_complete_count": fastqc_complete_count,
                    "fb_cpus_per_task": 1,
                    "fastqc_only_approx_sec": approx if approx is not None else "",
                    "fastqc_only_approx_basis": "whole fb job wall time; download and dumping skipped; FastQC completed"
                    if approx is not None
                    else "",
                }
            )
            enhanced_rows.append(enhanced)
    enhanced_path = FINAL / "timing/fb_get_srr_fastqc_log_summary_enhanced.tsv"
    if enhanced_rows:
        fieldnames = list(enhanced_rows[0].keys())
        with enhanced_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(enhanced_rows)
    out: dict[str, dict[str, object]] = {}
    for srr in set(wall) | set(dl) | set(log_paths):
        out[srr] = {
            "fb_get_srr_fastqc_wall_sec_min": min(wall[srr]) if wall[srr] else "",
            "fb_get_srr_fastqc_wall_sec_max": max(wall[srr]) if wall[srr] else "",
            "fb_get_srr_fastqc_wall_sec_median": median(wall[srr]) if wall[srr] else "",
            "fb_get_srr_fastqc_log_rows": len(log_paths[srr]),
            "fb_sra_https_download_sec_min": min(dl[srr]) if dl[srr] else "",
            "fb_sra_https_download_sec_max": max(dl[srr]) if dl[srr] else "",
            "fb_sra_https_download_sec_median": median(dl[srr]) if dl[srr] else "",
            "fb_fastqc_only_approx_sec_min": min(fastqc_approx[srr]) if fastqc_approx[srr] else "",
            "fb_fastqc_only_approx_sec_max": max(fastqc_approx[srr]) if fastqc_approx[srr] else "",
            "fb_fastqc_only_approx_sec_median": median(fastqc_approx[srr]) if fastqc_approx[srr] else "",
            "fb_fastqc_only_approx_log_rows": len(fastqc_approx[srr]),
            "fb_cpus_per_task": 1,
            "fb_log_path_example": log_paths[srr][0] if log_paths[srr] else "",
        }
    return out


def main() -> None:
    db = load_db()
    report_rows = load_manifest_report_rows()
    flat_timing = load_flattened_fastqc_timings()
    fb_timing = load_fb_timing_summary()

    fields = [
        "project_id",
        "srr_id",
        "report_source",
        "report_read_mate",
        "fastqc_data_path",
        "fastqc_zip_path",
        "db_status",
        "db_run_status",
        "db_fetch_lane",
        "db_size_gb",
        "db_sra_size_bytes",
        "db_fastq_size_bytes",
        "fastqc_version",
        "filename",
        "file_type",
        "encoding",
        "total_sequences",
        "total_bases",
        "sequences_flagged_poor_quality",
        "sequence_length",
        "gc_percent",
        "per_base_quality_status",
        "adapter_content_status",
        "sequence_duplication_status",
        "overrepresented_sequences_status",
        "mean_quality_all_bases",
        "mean_quality_first_bin",
        "mean_quality_terminal_bin",
        "mean_quality_last_10_bins",
        "median_quality_terminal_bin",
        "lowest_mean_quality_bin",
        "lowest_mean_quality_position",
        "fastqc_wall_sec_flattened_min",
        "fastqc_wall_sec_flattened_max",
        "fastqc_wall_sec_flattened_median",
        "fastqc_timing_rows",
        "fastqc_timing_core_assumption",
        "fastqc_timing_scope",
        "fb_get_srr_fastqc_wall_sec_min",
        "fb_get_srr_fastqc_wall_sec_max",
        "fb_get_srr_fastqc_wall_sec_median",
        "fb_get_srr_fastqc_log_rows",
        "fb_sra_https_download_sec_min",
        "fb_sra_https_download_sec_max",
        "fb_sra_https_download_sec_median",
        "fb_fastqc_only_approx_sec_min",
        "fb_fastqc_only_approx_sec_max",
        "fb_fastqc_only_approx_sec_median",
        "fb_fastqc_only_approx_log_rows",
        "fb_cpus_per_task",
        "fb_timing_core_note",
        "fb_timing_scope",
        "timing_preferred_wall_sec",
        "timing_preferred_scope",
        "timing_preferred_core_note",
        "fastqc_timing_path_example",
        "fb_log_path_example",
    ]

    rows: list[dict[str, object]] = []
    for srr in sorted(report_rows):
        manifest_row = report_rows[srr]
        data_path = Path(manifest_row["destination_path"])
        if not data_path.exists():
            continue
        db_row = db.get(srr, {})
        metrics = parse_fastqc(data_path)
        zip_path = str(data_path.parent) + ".zip"
        ftime = flat_timing.get(srr, {})
        btime = fb_timing.get(srr, {})
        preferred_wall = ""
        preferred_scope = ""
        preferred_core = ""
        if ftime.get("fastqc_wall_sec_flattened_max") not in (None, ""):
            preferred_wall = ftime.get("fastqc_wall_sec_flattened_max", "")
            preferred_scope = "FastQC stage only"
            preferred_core = "ja flattened timings; longest row used when multiple rows exist; treated as 1 core unless scheduler metadata says otherwise"
        elif btime.get("fb_fastqc_only_approx_sec_max") not in (None, ""):
            preferred_wall = btime.get("fb_fastqc_only_approx_sec_max", "")
            preferred_scope = "approximate fb FastQC-only wall time"
            preferred_core = "fb script requested --cpus-per-task=1; longest clean skip-download/skip-dump FastQC attempt used"
        elif btime.get("fb_get_srr_fastqc_wall_sec_max") not in (None, ""):
            preferred_wall = btime.get("fb_get_srr_fastqc_wall_sec_max", "")
            preferred_scope = "fb get-SRR plus conversion/download plus FastQC; not FastQC-only"
            preferred_core = "fb script requested --cpus-per-task=1; longest log used because shorter retry logs often represent failed attempts"

        row: dict[str, object] = {
            "project_id": manifest_row.get("project_id", db_row.get("project_id", "")),
            "srr_id": srr,
            "report_source": manifest_row.get("source", ""),
            "report_read_mate": read_mate_from_name(str(data_path)),
            "fastqc_data_path": str(data_path),
            "fastqc_zip_path": zip_path if Path(zip_path).exists() else "",
            "db_status": db_row.get("status", ""),
            "db_run_status": db_row.get("run_status", ""),
            "db_fetch_lane": db_row.get("fetch_lane", ""),
            "db_size_gb": db_row.get("size_gb", ""),
            "db_sra_size_bytes": db_row.get("sra_size_bytes", ""),
            "db_fastq_size_bytes": db_row.get("fastq_size_bytes", ""),
            **metrics,
            "fastqc_wall_sec_flattened_min": ftime.get("fastqc_wall_sec_flattened_min", ""),
            "fastqc_wall_sec_flattened_max": ftime.get("fastqc_wall_sec_flattened_max", ""),
            "fastqc_wall_sec_flattened_median": ftime.get("fastqc_wall_sec_flattened_median", ""),
            "fastqc_timing_rows": ftime.get("fastqc_timing_rows", ""),
            "fastqc_timing_core_assumption": "1 core for ja flattened timings; max duration used when multiple rows are present",
            "fastqc_timing_scope": "FastQC stage only",
            "fb_get_srr_fastqc_wall_sec_min": btime.get("fb_get_srr_fastqc_wall_sec_min", ""),
            "fb_get_srr_fastqc_wall_sec_max": btime.get("fb_get_srr_fastqc_wall_sec_max", ""),
            "fb_get_srr_fastqc_wall_sec_median": btime.get("fb_get_srr_fastqc_wall_sec_median", ""),
            "fb_get_srr_fastqc_log_rows": btime.get("fb_get_srr_fastqc_log_rows", ""),
            "fb_sra_https_download_sec_min": btime.get("fb_sra_https_download_sec_min", ""),
            "fb_sra_https_download_sec_max": btime.get("fb_sra_https_download_sec_max", ""),
            "fb_sra_https_download_sec_median": btime.get("fb_sra_https_download_sec_median", ""),
            "fb_fastqc_only_approx_sec_min": btime.get("fb_fastqc_only_approx_sec_min", ""),
            "fb_fastqc_only_approx_sec_max": btime.get("fb_fastqc_only_approx_sec_max", ""),
            "fb_fastqc_only_approx_sec_median": btime.get("fb_fastqc_only_approx_sec_median", ""),
            "fb_fastqc_only_approx_log_rows": btime.get("fb_fastqc_only_approx_log_rows", ""),
            "fb_cpus_per_task": btime.get("fb_cpus_per_task", ""),
            "fb_timing_core_note": "fb script requested --cpus-per-task=1; max wall time used across attempts; broad job timing is not FastQC-only",
            "fb_timing_scope": "whole get-SRR/FastQC job when available; includes download/conversion if performed in that attempt",
            "timing_preferred_wall_sec": preferred_wall,
            "timing_preferred_scope": preferred_scope,
            "timing_preferred_core_note": preferred_core,
            "fastqc_timing_path_example": ftime.get("fastqc_timing_path_example", ""),
            "fb_log_path_example": btime.get("fb_log_path_example", ""),
        }
        rows.append(row)

    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary_path = FINAL / "srr_fastqc_metrics_summary.tsv"
    summary = [
        ("rows", len(rows)),
        ("unique_srrs", len({r["srr_id"] for r in rows})),
        ("with_flattened_fastqc_stage_timing", sum(1 for r in rows if r["fastqc_wall_sec_flattened_median"] != "")),
        ("with_fb_fastqc_only_approx_timing", sum(1 for r in rows if r["fb_fastqc_only_approx_sec_max"] != "")),
        ("with_fb_get_srr_fastqc_wall_timing", sum(1 for r in rows if r["fb_get_srr_fastqc_wall_sec_median"] != "")),
        ("with_db_fastq_size_bytes", sum(1 for r in rows if r["db_fastq_size_bytes"] not in ("", None))),
        ("with_db_sra_size_bytes", sum(1 for r in rows if r["db_sra_size_bytes"] not in ("", None))),
    ]
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(summary)
    for metric, value in summary:
        print(f"{metric}\t{value}")
    print(f"wrote\t{OUT}")


if __name__ == "__main__":
    main()
