#!/usr/bin/env python3
"""Prune final_fastqc to one canonical FastQC report per DB SRR.

The srr_queue database is treated as the authority. Files not associated with
an SRR in the database are removed. For SRRs with multiple FastQC reports
(typically paired-end read 1/read 2), one canonical report is retained:
flattened source is preferred, then complete reports, then R1, then unpaired,
then R2.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import sqlite3
from collections import defaultdict
from pathlib import Path


BASE = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja")
OUT = BASE / "final_fastqc"
REPORTS = OUT / "reports"
MANIFEST = OUT / "manifest.tsv"
DB = Path("/pc2/users/o/omiks001/srr_queue.db")


def load_db_srrs() -> set[str]:
    con = sqlite3.connect(DB)
    try:
        return {row[0] for row in con.execute("select distinct srr_id from srr_queue")}
    finally:
        con.close()


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def srr_from_text(text: str) -> str:
    m = re.search(r"(SRR\d+)", text)
    return m.group(1) if m else ""


def report_base(path: Path) -> Path | None:
    """Return canonical extracted-report base for a FastQC report file."""
    name = path.name
    if name.endswith("_fastqc.zip"):
        return path.with_suffix("")
    if name.endswith("_fastqc.html"):
        return path.with_suffix("")
    for parent in [path, *path.parents]:
        if parent.name.endswith("_fastqc"):
            return parent
        if parent == REPORTS:
            break
    return None


def read_rank(base_name: str) -> int:
    if re.search(r"_1_fastqc$", base_name):
        return 3
    if re.search(r"_2_fastqc$", base_name):
        return 1
    return 2


def completeness_score(files: list[dict[str, str]]) -> int:
    dests = [Path(row["destination_path"]) for row in files]
    names = {p.name for p in dests}
    score = 0
    if any(p.name.endswith("_fastqc.zip") for p in dests):
        score += 4
    if "fastqc_data.txt" in names:
        score += 3
    if "summary.txt" in names:
        score += 2
    if "fastqc_report.html" in names or any(p.name.endswith("_fastqc.html") for p in dests):
        score += 1
    return score


def main() -> None:
    db_srrs = load_db_srrs()
    rows = read_manifest()

    backup = OUT / "manifest_before_db_prune.tsv"
    if not backup.exists():
        shutil.copy2(MANIFEST, backup)

    by_report: dict[Path, list[dict[str, str]]] = defaultdict(list)
    logs_or_other: list[dict[str, str]] = []
    for row in rows:
        dst = Path(row["destination_path"])
        base = report_base(dst)
        if base is None:
            logs_or_other.append(row)
        else:
            by_report[base].append(row)

    report_srr = {base: srr_from_text(str(base)) for base in by_report}

    report_groups_by_srr: dict[str, list[Path]] = defaultdict(list)
    for base, srr in report_srr.items():
        if srr in db_srrs:
            report_groups_by_srr[srr].append(base)

    chosen_by_srr: dict[str, Path] = {}
    for srr, bases in report_groups_by_srr.items():
        def key(base: Path) -> tuple[int, int, int, str]:
            group_rows = by_report[base]
            source_priority = max(int(r.get("source_priority") or 0) for r in group_rows)
            return (
                source_priority,
                completeness_score(group_rows),
                read_rank(base.name),
                str(base),
            )

        chosen_by_srr[srr] = max(bases, key=key)

    keep_files: set[Path] = set()
    kept_rows: list[dict[str, str]] = []
    prune_rows: list[dict[str, str]] = []

    for base, group_rows in by_report.items():
        srr = report_srr.get(base, "")
        keep = srr in chosen_by_srr and chosen_by_srr[srr] == base
        for row in group_rows:
            dst = Path(row["destination_path"])
            reason = "kept_canonical_db_srr_report" if keep else (
                "dropped_non_db_srr_report" if srr not in db_srrs else "dropped_extra_report_for_db_srr"
            )
            out_row = dict(row)
            out_row["prune_decision"] = "keep" if keep else "delete"
            out_row["prune_reason"] = reason
            prune_rows.append(out_row)
            if keep:
                keep_files.add(dst)
                kept_rows.append(row)

    for row in logs_or_other:
        dst = Path(row["destination_path"])
        srr = row.get("srr_id") or srr_from_text(str(dst))
        keep = srr in db_srrs
        out_row = dict(row)
        out_row["prune_decision"] = "keep" if keep else "delete"
        out_row["prune_reason"] = "kept_db_srr_log_or_auxiliary" if keep else "dropped_non_db_srr_log_or_auxiliary"
        prune_rows.append(out_row)
        if keep:
            keep_files.add(dst)
            kept_rows.append(row)

    deleted = 0
    for row in prune_rows:
        if row["prune_decision"] != "delete":
            continue
        dst = Path(row["destination_path"])
        if dst.exists() and dst.is_file():
            dst.unlink()
            deleted += 1

    # Remove empty directories below reports, deepest first.
    for root, dirs, files in os.walk(REPORTS, topdown=False):
        path = Path(root)
        if path != REPORTS:
            try:
                path.rmdir()
            except OSError:
                pass

    fieldnames = list(rows[0].keys()) + ["prune_decision", "prune_reason"]
    with (OUT / "db_prune_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(prune_rows)

    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept_rows)

    complete_reports = sum(1 for p in REPORTS.rglob("*_fastqc.zip"))
    represented_srrs = {srr_from_text(str(p)) for p in REPORTS.rglob("*_fastqc.zip")}
    represented_srrs.discard("")

    summary = [
        ("db_distinct_srrs", len(db_srrs)),
        ("db_srrs_with_fastqc_report", len(represented_srrs)),
        ("complete_fastqc_zip_reports_after_prune", complete_reports),
        ("manifest_rows_after_prune", len(kept_rows)),
        ("deleted_files", deleted),
        ("db_srrs_without_fastqc_report", len(db_srrs - represented_srrs)),
    ]
    with (OUT / "summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(summary)

    for metric, value in summary:
        print(f"{metric}\t{value}")


if __name__ == "__main__":
    main()
