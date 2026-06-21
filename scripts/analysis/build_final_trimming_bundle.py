#!/usr/bin/env python3
"""Build a consolidated trimming bundle from fb and ja outputs.

fb-derived trimming files are copied first. ja flattened/analysis trimming
files are copied second. Sequence/alignment payloads are excluded, and only SRRs
currently present in srr_queue.db are retained for per-SRR files.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import sqlite3
import subprocess
from pathlib import Path


BASE = Path("/pc2/users/o/omiks001/hpc-prf-omiks")
JA = BASE / "ja"
FB_RESULTS = BASE / "fb/omiks_project/results"
FLATTENED_STATS = JA / "flattened_trimmomatic_stats"
ANALYSIS_TECH = JA / "analysis/technical"
MULTIQC_DATA = JA / "multiqc_output/multiqc_data"
DB = Path("/pc2/users/o/omiks001/srr_queue.db")
OUT = JA / "final_trimming"

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


def is_payload(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in FORBIDDEN_EXTENSIONS)


def srr_from_text(text: str) -> str:
    match = re.search(r"(SRR\d+)", text)
    return match.group(1) if match else ""


def load_db_srrs() -> set[str]:
    con = sqlite3.connect(DB)
    try:
        return {row[0] for row in con.execute("select distinct srr_id from srr_queue")}
    finally:
        con.close()


def copy_file(src: Path, dst: Path) -> tuple[str, int]:
    if is_payload(src):
        raise ValueError(f"Refusing payload file: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    previous = dst.stat().st_size if dst.exists() else -1
    shutil.copy2(src, dst)
    return ("overwrote" if previous >= 0 else "created", previous)


def find_fb_files() -> list[Path]:
    if not FB_RESULTS.exists():
        return []
    cmd = [
        "find",
        str(FB_RESULTS),
        "-type",
        "f",
        "(",
        "-path",
        "*_trimmomatic*/*",
        "-o",
        "-path",
        "*/_out_err/run_trimming*",
        "-o",
        "-path",
        "*/_out_err/run_trimmomatic*",
        ")",
    ]
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
    return sorted(Path(line) for line in proc.stdout.splitlines() if line and not is_payload(Path(line)))


def find_ja_flattened_files() -> list[Path]:
    if not FLATTENED_STATS.exists():
        return []
    proc = subprocess.run(["find", str(FLATTENED_STATS), "-type", "f"], check=True, text=True, stdout=subprocess.PIPE)
    return sorted(Path(line) for line in proc.stdout.splitlines() if line and not is_payload(Path(line)))


def find_analysis_trimming_tables() -> list[Path]:
    files: list[Path] = []
    if ANALYSIS_TECH.exists():
        for path in ANALYSIS_TECH.glob("*.tsv"):
            name = path.name.lower()
            if "trimmomatic" in name or "retention" in name:
                files.append(path)
    if MULTIQC_DATA.exists():
        for name in ("trimmomatic_plot.txt", "multiqc_trimmomatic.txt"):
            path = MULTIQC_DATA / name
            if path.exists() and path.is_file():
                files.append(path)
    return sorted(files)


def rel_dest(src: Path, source: str) -> Path:
    if source == "fb_results":
        return OUT / "records" / "fb_results" / src.relative_to(FB_RESULTS)
    if source == "flattened_trimmomatic_stats":
        return OUT / "records" / "flattened_trimmomatic_stats" / src.relative_to(FLATTENED_STATS)
    if source == "analysis_technical":
        return OUT / "summary_tables" / src.name
    if source == "multiqc_data":
        return OUT / "summary_tables" / "multiqc_data" / src.name
    return OUT / "records" / source / src.name


def project_srr(src: Path, source: str) -> tuple[str, str]:
    if source == "fb_results":
        rel = src.relative_to(FB_RESULTS)
        return (rel.parts[0] if len(rel.parts) > 0 else "", rel.parts[1] if len(rel.parts) > 1 else srr_from_text(str(src)))
    if source == "flattened_trimmomatic_stats":
        rel = src.relative_to(FLATTENED_STATS)
        return (rel.parts[0] if len(rel.parts) > 0 else "", srr_from_text(str(src)))
    return ("", srr_from_text(str(src)))


def include_file(src: Path, source: str, db_srrs: set[str]) -> bool:
    if is_payload(src):
        return False
    srr = srr_from_text(str(src))
    # Per-SRR files must belong to the current DB. Aggregate summary tables do
    # not always contain one SRR in the path and are already pruned separately.
    if srr and srr not in db_srrs:
        return False
    return True


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    db_srrs = load_db_srrs()
    manifest: list[dict[str, object]] = []
    skipped_non_db = 0
    skipped_payload = 0

    source_batches = [
        ("fb_results", "trimming_log_or_small_output", find_fb_files()),
        ("flattened_trimmomatic_stats", "trimmomatic_stats", find_ja_flattened_files()),
        ("analysis_technical", "trimming_summary_table", [p for p in find_analysis_trimming_tables() if ANALYSIS_TECH in p.parents]),
        ("multiqc_data", "multiqc_trimming_table", [p for p in find_analysis_trimming_tables() if MULTIQC_DATA in p.parents]),
    ]

    for source, kind, files in source_batches:
        for src in files:
            if is_payload(src):
                skipped_payload += 1
                continue
            if not include_file(src, source, db_srrs):
                skipped_non_db += 1
                continue
            dst = rel_dest(src, source)
            action, previous = copy_file(src, dst)
            project, srr = project_srr(src, source)
            manifest.append(
                {
                    "source_priority": 1 if source == "fb_results" else 2,
                    "source": source,
                    "source_kind": kind,
                    "project_id": project,
                    "srr_id": srr,
                    "source_path": str(src),
                    "destination_path": str(dst),
                    "action": action,
                    "source_size_bytes": src.stat().st_size,
                    "previous_destination_size_bytes": previous if previous >= 0 else "",
                }
            )

    fields = [
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
    write_tsv(OUT / "manifest.tsv", manifest, fields)

    summary_rows = [
        {"metric": "db_srrs_current", "value": len(db_srrs)},
        {"metric": "manifest_rows", "value": len(manifest)},
        {"metric": "fb_files_copied", "value": sum(1 for r in manifest if r["source"] == "fb_results")},
        {"metric": "flattened_stats_files_copied", "value": sum(1 for r in manifest if r["source"] == "flattened_trimmomatic_stats")},
        {"metric": "analysis_summary_tables_copied", "value": sum(1 for r in manifest if r["source"] == "analysis_technical")},
        {"metric": "multiqc_tables_copied", "value": sum(1 for r in manifest if r["source"] == "multiqc_data")},
        {"metric": "overwrites_by_later_layers", "value": sum(1 for r in manifest if r["action"] == "overwrote")},
        {"metric": "skipped_payload_files", "value": skipped_payload},
        {"metric": "skipped_non_db_srr_files", "value": skipped_non_db},
        {"metric": "unique_srrs_in_manifest", "value": len({r["srr_id"] for r in manifest if r["srr_id"]})},
    ]
    write_tsv(OUT / "summary.tsv", summary_rows, ["metric", "value"])
    for row in summary_rows:
        print(f"{row['metric']}\t{row['value']}")


if __name__ == "__main__":
    main()
