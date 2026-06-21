#!/usr/bin/env python3
"""
Extract per-SRR QC features from flattened FastQC directories.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path


FASTQC_BASE = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/flattened_fastqc_raw")
OUT_DIR = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis")

HEADER = [
    "project_id",
    "SRR_ID",
    "sequence_depth",
    "Q_mean",
    "Q_median",
    "Q_max",
    "Q_min",
    "read_length_mean",
    "read_length_median",
    "read_length_max",
    "read_length_min",
    "frac_below_q20",
    "frac_below_q30",
    "tail_quality_decay",
    "adapter_rate",
    "duplication_rate",
    "n_content",
    "gc_content",
    "gc_deviation",
]


def parse_range_width(token: str) -> int:
    token = token.strip()
    if "-" in token:
        a, b = token.split("-", 1)
        return int(b) - int(a) + 1
    return 1


def parse_range_midpoint(token: str) -> float:
    token = token.strip()
    if "-" in token:
        a, b = token.split("-", 1)
        return (int(a) + int(b)) / 2
    return float(token)


def parse_fastqc_data(path: Path) -> dict[str, float | str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    basic = {}
    sections: dict[str, list[list[str]]] = defaultdict(list)
    current = None
    for line in lines:
        if line.startswith(">>"):
            if line.startswith(">>END_MODULE"):
                current = None
                continue
            current = line[2:].split("\t", 1)[0]
            continue
        if current == "Basic Statistics":
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                basic[parts[0]] = parts[1]
        elif current:
            if line.startswith("#") or not line.strip():
                continue
            sections[current].append(line.split("\t"))

    total_sequences = float(basic.get("Total Sequences", "0") or 0)
    gc_content = float(basic.get("%GC", "0") or 0)

    # Per-base quality features
    quality_rows = sections.get("Per base sequence quality", [])
    q_weight = 0
    q_mean_sum = q_median_sum = 0.0
    q_max = float("-inf")
    q_min = float("inf")
    q20_weight = q30_weight = 0
    quality_profile = []
    for row in quality_rows:
        width = parse_range_width(row[0])
        mean = float(row[1])
        median = float(row[2])
        low_q = float(row[3])
        high_q = float(row[4])
        q_weight += width
        q_mean_sum += mean * width
        q_median_sum += median * width
        q_max = max(q_max, high_q)
        q_min = min(q_min, low_q)
        if mean < 20:
            q20_weight += width
        if mean < 30:
            q30_weight += width
        quality_profile.extend([mean] * width)

    q_mean = q_mean_sum / q_weight if q_weight else None
    q_median = q_median_sum / q_weight if q_weight else None
    frac_below_q20 = q20_weight / q_weight if q_weight else None
    frac_below_q30 = q30_weight / q_weight if q_weight else None
    if quality_profile:
        head = quality_profile[: min(10, len(quality_profile))]
        tail = quality_profile[-min(10, len(quality_profile)) :]
        tail_quality_decay = (sum(head) / len(head)) - (sum(tail) / len(tail))
    else:
        tail_quality_decay = None

    # Read length distribution
    length_rows = sections.get("Sequence Length Distribution", [])
    total_len_count = 0
    weighted_len_sum = 0.0
    lengths = []
    length_min = float("inf")
    length_max = float("-inf")
    for row in length_rows:
        length = parse_range_midpoint(row[0])
        count = float(row[1])
        total_len_count += count
        weighted_len_sum += length * count
        lengths.append((length, count))
        length_min = min(length_min, length)
        length_max = max(length_max, length)
    if lengths and total_len_count:
        read_length_mean = weighted_len_sum / total_len_count
        half = total_len_count / 2
        cum = 0.0
        read_length_median = read_length_mean
        for length, count in sorted(lengths):
            cum += count
            if cum >= half:
                read_length_median = length
                break
    else:
        read_length_mean = read_length_median = None
        length_min = length_max = None

    # Adapter content
    adapter_rows = sections.get("Adapter Content", [])
    adapter_rate = 0.0
    for row in adapter_rows:
        vals = [float(v) for v in row[1:] if v not in {"", "NA"}]
        if vals:
            adapter_rate = max(adapter_rate, max(vals))

    # Duplication
    dup_rows = sections.get("Sequence Duplication Levels", [])
    duplication_rate = None
    for row in dup_rows:
        if row[0] == "#Total Deduplicated Percentage" and len(row) > 1:
            duplication_rate = 100.0 - float(row[1])
            break
    if duplication_rate is None:
        for raw in lines:
            if raw.startswith("#Total Deduplicated Percentage"):
                parts = raw.split("\t")
                if len(parts) > 1:
                    duplication_rate = 100.0 - float(parts[1])
                break

    # N content
    n_rows = sections.get("Per base N content", [])
    n_weight = 0
    n_sum = 0.0
    for row in n_rows:
        width = parse_range_width(row[0])
        val = float(row[1])
        n_weight += width
        n_sum += val * width
    n_content = n_sum / n_weight if n_weight else 0.0

    return {
        "sequence_depth": total_sequences,
        "Q_mean": q_mean,
        "Q_median": q_median,
        "Q_max": q_max if q_weight else None,
        "Q_min": q_min if q_weight else None,
        "read_length_mean": read_length_mean,
        "read_length_median": read_length_median,
        "read_length_max": length_max,
        "read_length_min": length_min,
        "frac_below_q20": frac_below_q20,
        "frac_below_q30": frac_below_q30,
        "tail_quality_decay": tail_quality_decay,
        "adapter_rate": adapter_rate,
        "duplication_rate": duplication_rate,
        "n_content": n_content,
        "gc_content": gc_content,
        "gc_deviation": abs(gc_content - 50.0),
    }


def weighted_merge(rows: list[dict[str, float | str]]) -> dict[str, str]:
    total_depth = sum(float(r["sequence_depth"]) for r in rows if r["sequence_depth"] is not None)
    out = {
        "sequence_depth": str(int(round(total_depth))) if total_depth else "NA",
    }
    weighted_fields = [
        "Q_mean",
        "Q_median",
        "read_length_mean",
        "read_length_median",
        "frac_below_q20",
        "frac_below_q30",
        "tail_quality_decay",
        "adapter_rate",
        "duplication_rate",
        "n_content",
        "gc_content",
        "gc_deviation",
    ]
    max_fields = ["Q_max", "read_length_max"]
    min_fields = ["Q_min", "read_length_min"]

    for field in weighted_fields:
        vals = [(float(r["sequence_depth"]), float(r[field])) for r in rows if r.get(field) is not None]
        if vals and total_depth:
            out[field] = f"{sum(w * v for w, v in vals) / total_depth:.4f}"
        else:
            out[field] = "NA"

    for field in max_fields:
        vals = [float(r[field]) for r in rows if r.get(field) is not None]
        out[field] = f"{max(vals):.4f}" if vals else "NA"
    for field in min_fields:
        vals = [float(r[field]) for r in rows if r.get(field) is not None]
        out[field] = f"{min(vals):.4f}" if vals else "NA"

    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_filter", nargs="?", default="")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_file = OUT_DIR / "per_srr_quality.tsv"

    grouped: dict[tuple[str, str], list[dict[str, float | str]]] = defaultdict(list)
    projects = sorted(p for p in FASTQC_BASE.glob("*/") if p.is_dir())
    for project_dir in projects:
        project = project_dir.name
        if args.project_filter and project != args.project_filter:
            continue
        print(f"Processing {project}...")
        for fastqc_dir in sorted(project_dir.glob("*_fastqc")):
            data_file = fastqc_dir / "fastqc_data.txt"
            if not data_file.exists():
                continue
            srr = fastqc_dir.name.removesuffix("_fastqc")
            if srr.endswith("_1") or srr.endswith("_2"):
                srr = srr[:-2]
            grouped[(project, srr)].append(parse_fastqc_data(data_file))

    with out_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        for (project, srr) in sorted(grouped):
            merged = weighted_merge(grouped[(project, srr)])
            row = {"project_id": project, "SRR_ID": srr}
            row.update(merged)
            writer.writerow(row)

    print(f"Done. Total FastQC entries: {len(grouped)}, Unique SRRs: {len(grouped)}")
    print(f"Output: {out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
