#!/usr/bin/env python3
"""Archive FastQC files pruned from final_fastqc.

The pruning step removed non-canonical/non-DB copies from final_fastqc. The
original source files were not deleted. This script reconstructs an archive of
the pruned files from source_path values recorded in db_prune_manifest.tsv.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


FINAL = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_fastqc")
PRUNE_MANIFEST = FINAL / "db_prune_manifest.tsv"
ARCHIVE = FINAL / "archive_pruned_reports"

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


def forbidden(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in FORBIDDEN_EXTENSIONS)


def archive_destination(source: Path, row: dict[str, str]) -> Path:
    source_name = row.get("source", "unknown_source") or "unknown_source"
    project = row.get("project_id", "unknown_project") or "unknown_project"
    srr = row.get("srr_id", "unknown_srr") or "unknown_srr"
    try:
        if "hpc-prf-omiks/fb/omiks_project/results" in str(source):
            marker = "hpc-prf-omiks/fb/omiks_project/results"
            rel = Path(str(source).split(marker, 1)[1].lstrip("/"))
        elif "hpc-prf-omiks/ja/flattened_fastqc_raw" in str(source):
            marker = "hpc-prf-omiks/ja/flattened_fastqc_raw"
            rel = Path(str(source).split(marker, 1)[1].lstrip("/"))
        else:
            rel = Path(project) / srr / source.name
    except Exception:
        rel = Path(project) / srr / source.name
    return ARCHIVE / source_name / rel


def main() -> None:
    if not PRUNE_MANIFEST.exists():
        raise SystemExit(f"Missing prune manifest: {PRUNE_MANIFEST}")

    rows_out: list[dict[str, str]] = []
    copied = 0
    missing_source = 0
    skipped_payload = 0
    already_present = 0

    with PRUNE_MANIFEST.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("prune_decision") != "delete":
                continue
            src = Path(row.get("source_path", ""))
            out_row = dict(row)
            if not src.exists():
                missing_source += 1
                out_row["archive_status"] = "missing_source"
                out_row["archive_path"] = ""
                rows_out.append(out_row)
                continue
            if forbidden(src):
                skipped_payload += 1
                out_row["archive_status"] = "skipped_payload"
                out_row["archive_path"] = ""
                rows_out.append(out_row)
                continue
            dst = archive_destination(src, row)
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists():
                already_present += 1
                status = "already_present"
            else:
                shutil.copy2(src, dst)
                copied += 1
                status = "copied"
            out_row["archive_status"] = status
            out_row["archive_path"] = str(dst)
            rows_out.append(out_row)

    manifest_path = FINAL / "archive_pruned_manifest.tsv"
    fieldnames = list(rows_out[0].keys()) if rows_out else ["archive_status", "archive_path"]
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)

    summary = [
        ("archive_dir", str(ARCHIVE)),
        ("deleted_manifest_rows_considered", len(rows_out)),
        ("copied", copied),
        ("already_present", already_present),
        ("missing_source", missing_source),
        ("skipped_payload", skipped_payload),
    ]
    with (FINAL / "archive_pruned_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(summary)
    for metric, value in summary:
        print(f"{metric}\t{value}")


if __name__ == "__main__":
    main()
