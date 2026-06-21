#!/usr/bin/env python3
"""Archive non-Illumina and suspicious SRRs from on-disk records and DB.

This script moves per-SRR result records out of active result directories into
an archive and deletes the same SRR accessions from srr_queue.db. It does not
remove raw source data outside the ja result tree, and it keeps all cleaning
strategy outputs for retained SRRs.
"""

from __future__ import annotations

import csv
import os
import shutil
import sqlite3
from datetime import datetime
from pathlib import Path

import pandas as pd


JA = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja")
ANALYSIS = JA / "analysis"
TECH = ANALYSIS / "technical"
DB = Path("/pc2/users/o/omiks001/srr_queue.db")

SUSPICIOUS_SRRS = set(
    """
SRR35798687 SRR35798688 SRR35798695 SRR35798696 SRR35798697 SRR35798698
SRR6334449 SRR7058572 SRR7058573 SRR7058574 SRR7058575 SRR7058576
SRR7058577 SRR7058578 SRR7058579 SRR7058580 SRR19842866 SRR17165227
""".split()
)


def load_targets() -> pd.DataFrame:
    platform = pd.read_csv(TECH / "srr_platform_check.tsv", sep="\t")
    non_illumina = platform.loc[platform["is_illumina"].ne("yes"), [
        "run_accession",
        "study_accession",
        "instrument_platform",
        "library_strategy",
    ]].copy()
    non_illumina = non_illumina.rename(columns={"run_accession": "srr_id", "study_accession": "project_id"})
    non_illumina["prune_reason"] = "non_illumina"

    suspicious = pd.DataFrame({"srr_id": sorted(SUSPICIOUS_SRRS)})
    suspicious["project_id"] = ""
    suspicious["instrument_platform"] = ""
    suspicious["library_strategy"] = ""
    suspicious["prune_reason"] = "suspicious"

    targets = pd.concat([non_illumina, suspicious], ignore_index=True)
    targets = (
        targets.groupby("srr_id", as_index=False)
        .agg(
            project_id=("project_id", lambda x: ";".join(sorted(set(v for v in x if v)))),
            instrument_platform=("instrument_platform", lambda x: ";".join(sorted(set(v for v in x if v)))),
            library_strategy=("library_strategy", lambda x: ";".join(sorted(set(v for v in x if v)))),
            prune_reason=("prune_reason", lambda x: ";".join(sorted(set(x)))),
        )
        .sort_values("srr_id")
    )
    return targets


def db_rows_for_targets(target_srrs: set[str]) -> list[dict[str, object]]:
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    try:
        placeholders = ",".join("?" for _ in target_srrs)
        if not placeholders:
            return []
        rows = con.execute(f"select * from srr_queue where srr_id in ({placeholders})", sorted(target_srrs)).fetchall()
        return [dict(r) for r in rows]
    finally:
        con.close()


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def backup_db(archive_dir: Path) -> Path:
    backup_path = archive_dir / f"srr_queue.before_filtered_srr_prune.bak"
    shutil.copy2(DB, backup_path)
    return backup_path


def top_level_scan_roots() -> list[Path]:
    roots = []
    for path in JA.iterdir():
        if not path.is_dir():
            continue
        name = path.name
        if name.startswith("flattened_") or name in {"final_fastqc", "multiqc_output"}:
            roots.append(path)
    return sorted(roots)


def path_contains_srr(path: Path, target_srrs: set[str]) -> str | None:
    text_parts = path.parts
    for srr in target_srrs:
        if any(srr in part for part in text_parts):
            return srr
    return None


def relative_archive_destination(src: Path, archive_dir: Path) -> Path:
    rel = src.relative_to(JA)
    return archive_dir / "disk_records" / rel


def move_matching_records(target_srrs: set[str], archive_dir: Path) -> list[dict[str, object]]:
    move_rows: list[dict[str, object]] = []
    for root in top_level_scan_roots():
        if not root.exists():
            continue
        # Walk top-down so a matched SRR directory is moved as one unit.
        for current, dirs, files in os.walk(root, topdown=True):
            current_path = Path(current)
            matched_dir_srr = path_contains_srr(current_path, target_srrs)
            if matched_dir_srr and current_path != root:
                dst = relative_archive_destination(current_path, archive_dir)
                dst.parent.mkdir(parents=True, exist_ok=True)
                if dst.exists():
                    suffix = datetime.now().strftime("%H%M%S%f")
                    dst = dst.with_name(dst.name + f".dup_{suffix}")
                shutil.move(str(current_path), str(dst))
                move_rows.append({
                    "srr_id": matched_dir_srr,
                    "record_type": "directory",
                    "source_path": str(current_path),
                    "archive_path": str(dst),
                })
                dirs[:] = []
                continue

            kept_dirs = []
            for d in dirs:
                dpath = current_path / d
                srr = path_contains_srr(dpath, target_srrs)
                if srr:
                    dst = relative_archive_destination(dpath, archive_dir)
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    if dst.exists():
                        suffix = datetime.now().strftime("%H%M%S%f")
                        dst = dst.with_name(dst.name + f".dup_{suffix}")
                    shutil.move(str(dpath), str(dst))
                    move_rows.append({
                        "srr_id": srr,
                        "record_type": "directory",
                        "source_path": str(dpath),
                        "archive_path": str(dst),
                    })
                else:
                    kept_dirs.append(d)
            dirs[:] = kept_dirs

            for f in files:
                fpath = current_path / f
                srr = path_contains_srr(fpath, target_srrs)
                if not srr:
                    continue
                dst = relative_archive_destination(fpath, archive_dir)
                dst.parent.mkdir(parents=True, exist_ok=True)
                if dst.exists():
                    suffix = datetime.now().strftime("%H%M%S%f")
                    dst = dst.with_name(dst.name + f".dup_{suffix}")
                shutil.move(str(fpath), str(dst))
                move_rows.append({
                    "srr_id": srr,
                    "record_type": "file",
                    "source_path": str(fpath),
                    "archive_path": str(dst),
                })
    return move_rows


def delete_db_rows(target_srrs: set[str]) -> int:
    con = sqlite3.connect(DB)
    try:
        before = con.execute("select count(*) from srr_queue").fetchone()[0]
        con.executemany("delete from srr_queue where srr_id = ?", [(srr,) for srr in sorted(target_srrs)])
        con.commit()
        after = con.execute("select count(*) from srr_queue").fetchone()[0]
        return before - after
    finally:
        con.close()


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_dir = JA / "archive" / f"filtered_srrs_{timestamp}"
    archive_dir.mkdir(parents=True, exist_ok=True)

    targets = load_targets()
    target_srrs = set(targets["srr_id"])
    targets.to_csv(archive_dir / "target_srrs.tsv", sep="\t", index=False)

    db_rows = db_rows_for_targets(target_srrs)
    write_tsv(archive_dir / "deleted_db_rows.tsv", db_rows)
    db_target_srrs = {row["srr_id"] for row in db_rows}

    backup_path = backup_db(archive_dir)
    move_rows = move_matching_records(target_srrs, archive_dir)
    write_tsv(archive_dir / "move_manifest.tsv", move_rows)
    deleted_db_rows = delete_db_rows(target_srrs)

    summary = [
        {"metric": "archive_dir", "value": str(archive_dir)},
        {"metric": "db_backup", "value": str(backup_path)},
        {"metric": "target_srrs_total", "value": len(target_srrs)},
        {"metric": "target_srrs_in_db", "value": len(db_target_srrs)},
        {"metric": "target_srrs_non_illumina", "value": int(targets["prune_reason"].str.contains("non_illumina").sum())},
        {"metric": "target_srrs_suspicious", "value": int(targets["prune_reason"].str.contains("suspicious").sum())},
        {"metric": "disk_records_moved", "value": len(move_rows)},
        {"metric": "db_rows_deleted", "value": deleted_db_rows},
    ]
    write_tsv(archive_dir / "summary.tsv", summary, ["metric", "value"])

    for row in summary:
        print(f"{row['metric']}\t{row['value']}")


if __name__ == "__main__":
    main()
