#!/usr/bin/env python3
"""Remove archived/filtered SRR rows from active TSV tables.

Per-SRR files are moved by archive_and_prune_filtered_srrs.py. This companion
script removes rows for the same SRRs from active aggregate TSV tables and
stores both backups and removed rows in the same archive.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

import pandas as pd


JA = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja")
ARCHIVE_ROOT = JA / "archive"

SRR_COLUMNS = [
    "srr_id",
    "SRR_ID",
    "run_accession",
    "Run",
    "sample",
    "sample_id",
]


def latest_archive() -> Path:
    archives = sorted(ARCHIVE_ROOT.glob("filtered_srrs_*"))
    if not archives:
        raise SystemExit("No filtered_srrs_* archive found")
    return archives[-1]


def load_targets(archive: Path) -> set[str]:
    targets = pd.read_csv(archive / "target_srrs.tsv", sep="\t")
    return set(targets["srr_id"].astype(str))


def candidate_tsvs() -> list[Path]:
    roots = [
        JA / "final_fastqc",
        JA / "analysis" / "technical",
        JA / "analysis" / "deseq2_metadata",
    ]
    files: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*.tsv"):
            if "/archive/" in str(path):
                continue
            if "archive_pruned_reports" in str(path) or "archive_filtered" in str(path):
                continue
            files.append(path)
    return sorted(files)


def row_matches(row: dict[str, str], target_srrs: set[str]) -> bool:
    for col in SRR_COLUMNS:
        val = row.get(col)
        if val and val in target_srrs:
            return True
    # Some manifest/source-path tables contain SRRs in paths rather than a
    # normalized accession column.
    for val in row.values():
        if not val:
            continue
        for srr in target_srrs:
            if srr in val:
                return True
    return False


def prune_tsv(path: Path, target_srrs: set[str], archive: Path) -> dict[str, object] | None:
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = reader.fieldnames
            if not fieldnames:
                return None
            kept = []
            removed = []
            for row in reader:
                (removed if row_matches(row, target_srrs) else kept).append(row)
    except UnicodeDecodeError:
        return None

    if not removed:
        return None

    rel = path.relative_to(JA)
    backup = archive / "table_backups" / rel
    removed_path = archive / "removed_table_rows" / rel
    backup.parent.mkdir(parents=True, exist_ok=True)
    removed_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup)

    with removed_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(removed)

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept)

    return {
        "table_path": str(path),
        "backup_path": str(backup),
        "removed_rows_path": str(removed_path),
        "kept_rows": len(kept),
        "removed_rows": len(removed),
    }


def main() -> None:
    archive = latest_archive()
    target_srrs = load_targets(archive)
    summaries = []
    for path in candidate_tsvs():
        result = prune_tsv(path, target_srrs, archive)
        if result:
            summaries.append(result)

    summary_path = archive / "table_prune_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        fieldnames = ["table_path", "backup_path", "removed_rows_path", "kept_rows", "removed_rows"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)

    print(f"archive\t{archive}")
    print(f"tables_pruned\t{len(summaries)}")
    print(f"rows_removed\t{sum(int(r['removed_rows']) for r in summaries)}")
    print(f"summary\t{summary_path}")


if __name__ == "__main__":
    main()
