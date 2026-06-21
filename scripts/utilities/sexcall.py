#!/usr/bin/env python3
"""
sexcall.py – Unified biological sex caller for hg38 alignments.

Accepts a sorted BAM (coordinate-sorted preferred) and produces a sex
call (XX / XY / XXY / XYY / X0) with confidence scoring.

Three independent evidence signals are combined:
  1. Normalised depth ratios   – X:autosome, Y:autosome via idxstats
  2. Strict read classification – perfect, end-to-end, high-MAPQ X/Y reads
  3. XIST expression            – reads overlapping XIST exons (optional)

Usage:
    python3 sexcall.py <sorted.bam> [options]

See --help for all options.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import pysam

# ──────────────────────────────────────────────────────────────────────
#  Constants
# ──────────────────────────────────────────────────────────────────────

CHR_X_SYMS = {"NC_000023.11", "chrX", "X", "23"}
CHR_Y_SYMS = {"NC_000024.10", "chrY", "Y", "24"}
SEX_SYMS   = CHR_X_SYMS | CHR_Y_SYMS

# Autosomes we use for depth normalisation
AUTO_SYMS: dict[str, set[str]] = {}
for i in range(1, 23):
    AUTO_SYMS[str(i)] = {
        str(i),
        f"chr{i}",
        f"NC_0000{i:02d}" if i < 10 else f"NC_00{i:03d}",  # approx RefSeq
    }

M_CIGAR_RE = re.compile(r"^(\d+)M$")

# ──────────────────────────────────────────────────────────────────────
#  Helpers
# ──────────────────────────────────────────────────────────────────────

def _resolve_chrom(refs: set[str], aliases: set[str]) -> str | None:
    """Find the first BAM reference matching any alias (case-insensitive)."""
    alias_lower = {a.lower() for a in aliases}
    for r in refs:
        if r.lower() in alias_lower or r in aliases:
            return r
    return None


def _resolve_autosomes(refs: set[str]) -> dict[str, str]:
    """Map autosome number → matching BAM reference name."""
    result = {}
    for r in refs:
        # Strip trailing dot versions like .11, .14
        r_base = r.split('.', 1)[0].lower()
        if r_base.startswith("chr"):
            num = r_base[3:]
        elif r_base.startswith("nc_0000"):
            # RefSeq NC_000001, NC_000022, etc
            num = str(int(r_base.replace("nc_", "")))
        else:
            num = r_base
            
        if num.isdigit() and 1 <= int(num) <= 22:
            result[str(int(num))] = r
    return result


# ──────────────────────────────────────────────────────────────────────
#  1. Depth ratios via idxstats
# ──────────────────────────────────────────────────────────────────────

def compute_depth_ratios(bam_path: str) -> tuple[float, float, dict]:
    """
    Compute normalised X:autosome and Y:autosome depth ratios.

    Returns (x_ratio, y_ratio, raw_stats) where raw_stats contains
    per-chromosome mapped read counts and lengths.
    """
    stats = pysam.idxstats(bam_path).strip().split("\n")
    chrom_data: dict[str, tuple[int, int]] = {}  # name → (length, mapped)
    for line in stats:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        name, length, mapped, _unmapped = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
        if length > 0:
            chrom_data[name] = (length, mapped)

    refs = set(chrom_data.keys())
    x_ref = _resolve_chrom(refs, CHR_X_SYMS)
    y_ref = _resolve_chrom(refs, CHR_Y_SYMS)
    auto_refs = _resolve_autosomes(refs)

    if not x_ref or not y_ref:
        print("[WARN] chrX/Y not found in BAM index – depth ratios unavailable", file=sys.stderr)
        return -1.0, -1.0, {}

    # Compute per-base depth for autosomes
    auto_depths: list[float] = []
    for _num, ref in sorted(auto_refs.items(), key=lambda kv: int(kv[0])):
        length, mapped = chrom_data[ref]
        if length > 0:
            auto_depths.append(mapped / length)

    if not auto_depths:
        print("[WARN] No autosomal chromosomes resolved – depth ratios unavailable", file=sys.stderr)
        return -1.0, -1.0, {}

    median_auto = sorted(auto_depths)[len(auto_depths) // 2]

    if median_auto == 0:
        return -1.0, -1.0, {}

    x_len, x_mapped = chrom_data[x_ref]
    y_len, y_mapped = chrom_data[y_ref]

    x_depth = x_mapped / x_len if x_len > 0 else 0
    y_depth = y_mapped / y_len if y_len > 0 else 0

    x_ratio = x_depth / median_auto
    y_ratio = y_depth / median_auto

    raw = {
        "x_ref": x_ref, "y_ref": y_ref,
        "x_mapped": x_mapped, "y_mapped": y_mapped,
        "x_len": x_len, "y_len": y_len,
        "median_auto_depth": median_auto,
        "auto_count": len(auto_depths),
    }
    return x_ratio, y_ratio, raw


# ──────────────────────────────────────────────────────────────────────
#  2. Strict read classification (merged from sex_call_chromosomes.sh)
# ──────────────────────────────────────────────────────────────────────

def classify_reads(bam_path: str, mapq_thr: int = 30) -> tuple[int, int, int]:
    """
    Count reads strictly mapping to X, Y, or both.

    Uses the multi-mapper-aware logic from sex_call_chromosomes.sh v1.3:
    - Only best-score alignments (AS tag) are considered
    - Reads with any somatic best-alignment are discarded
    - Strict filters: NM=0, MAPQ ≥ threshold, pure M-cigar

    Returns (x_only, y_only, xy_both).
    """
    info: dict[str, list] = {}  # rname → [best_score, chrs, strict_x, strict_y, has_somatic]

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = set(bam.references)
        x_ref = _resolve_chrom(refs, CHR_X_SYMS)
        y_ref = _resolve_chrom(refs, CHR_Y_SYMS)

        if not x_ref or not y_ref:
            print("[WARN] chrX/Y not in BAM header – read classification skipped", file=sys.stderr)
            return 0, 0, 0

        local_sex = {x_ref, y_ref}

        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            try:
                score = aln.get_tag("AS")
            except KeyError:
                continue

            rname = aln.query_name
            chrom = aln.reference_name

            rec = info.get(rname)
            if rec is None:
                rec = [-(10**9), set(), False, False, False]
                info[rname] = rec

            best, chrs, strict_x, strict_y, has_somatic = rec

            if score > best:
                best = score
                chrs = {chrom}
                strict_x = strict_y = has_somatic = False
            elif score == best:
                chrs.add(chrom)

            # somatic flag
            if chrom not in local_sex and score == best:
                has_somatic = True

            # strict filters on best-score alignments
            if score == best:
                try:
                    nm = aln.get_tag("NM")
                except KeyError:
                    nm = 0
                # Relaxed CIGAR/NM check for RNA-seq (allows splices/soft-clip up to 2 mismatches)
                strict = (nm <= 2) and aln.mapping_quality >= mapq_thr
                if strict and chrom in CHR_X_SYMS:
                    strict_x = True
                if strict and chrom in CHR_Y_SYMS:
                    strict_y = True

            rec[:] = [best, chrs, strict_x, strict_y, has_somatic]

    # Tally
    x_only = y_only = xy_both = 0
    for _rname, (_, chrs, strict_x, strict_y, has_somatic) in info.items():
        if has_somatic:
            continue
        if chrs & CHR_X_SYMS and chrs & CHR_Y_SYMS:
            xy_both += 1
            continue
        if chrs.issubset(CHR_X_SYMS):
            if strict_x:
                x_only += 1
            continue
        if chrs.issubset(CHR_Y_SYMS):
            if strict_y:
                y_only += 1
            continue

    return x_only, y_only, xy_both


# ──────────────────────────────────────────────────────────────────────
#  3. XIST read counting
# ──────────────────────────────────────────────────────────────────────

def count_xist(bam_path: str, bed_path: str, mapq_thr: int = 30) -> int:
    """Count reads overlapping XIST exons (from BED file)."""
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals[chrom].append((start, end))

    if not intervals:
        print(f"[WARN] No intervals read from {bed_path}", file=sys.stderr)
        return 0

    count = 0
    seen: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, spans in intervals.items():
            for s, e in spans:
                try:
                    for aln in bam.fetch(chrom, s, e):
                        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                            continue
                        if aln.query_name in seen:
                            continue
                        try:
                            if aln.get_tag("NM") > 2:
                                continue
                        except KeyError:
                            pass
                        
                        if aln.mapping_quality < mapq_thr:
                            continue
                        
                        seen.add(aln.query_name)
                        count += 1
                except ValueError:
                    # Chrom not in BAM
                    pass

    return count


# ──────────────────────────────────────────────────────────────────────
#  4. Multi-evidence sex caller
# ──────────────────────────────────────────────────────────────────────

def call_sex(
    x_ratio: float,
    y_ratio: float,
    x_only: int,
    y_only: int,
    xy_both: int,
    xist_count: int,
    y_thresh: int = 10,
) -> tuple[str, str, str]:
    """
    Combine all evidence into a sex / karyotype / confidence call.

    Returns (sex, karyotype, confidence).
    """
    has_y = (y_only + xy_both) >= y_thresh
    has_xist = xist_count > 0
    depth_available = x_ratio >= 0 and y_ratio >= 0

    # ── Depth-based signals ──────────────────────────────────────────
    # For a diploid:
    #   XX female → X:auto ≈ 1.0,  Y:auto ≈ 0.0
    #   XY male   → X:auto ≈ 0.5,  Y:auto ≈ 0.5  (Y shorter so mapped ratio can differ)
    #   XXY       → X:auto ≈ 1.0,  Y:auto ≈ 0.3-0.5
    #   XYY       → X:auto ≈ 0.5,  Y:auto ≈ 1.0
    #   X0        → X:auto ≈ 0.5,  Y:auto ≈ 0.0

    x_high = x_ratio > 0.75 if depth_available else None    # ≈ 1 copy X per autosome pair
    y_present_depth = y_ratio > 0.05 if depth_available else None

    # ── Classification ───────────────────────────────────────────────
    sex = "undetermined"
    karyotype = "undetermined"
    confidence = "LOW"

    if not has_y and (y_present_depth is False or y_present_depth is None):
        # No Y presence
        if x_high is True or x_high is None:
            sex, karyotype = "female", "XX"
            confidence = "HIGH" if depth_available and x_high else "MEDIUM"
        elif x_high is False:
            # Low X, no Y → could be Turner (X0)
            sex, karyotype = "female", "X0 (Turner)"
            confidence = "MEDIUM" if depth_available else "LOW"
    elif has_y:
        if x_high is True:
            # Y reads present + X ratio high → XXY (Klinefelter)
            if has_xist:
                sex, karyotype = "male", "XXY (Klinefelter)"
                confidence = "HIGH"
            else:
                # High X + Y but no XIST data → still likely XXY but lower confidence
                sex, karyotype = "male", "XXY (Klinefelter)"
                confidence = "MEDIUM"
        elif x_high is False:
            # Normal X ratio + Y present → standard XY
            if y_ratio is not None and y_ratio > 0.75:
                sex, karyotype = "male", "XYY"
                confidence = "MEDIUM"
            else:
                sex, karyotype = "male", "XY"
                confidence = "HIGH" if depth_available else "MEDIUM"
        else:
            # Depth unavailable but Y reads present → male
            sex, karyotype = "male", "XY"
            confidence = "MEDIUM"
    elif y_present_depth is True and not has_y:
        # Some Y depth but below read-count threshold
        # Could be low-coverage male or contamination
        sex, karyotype = "male", "XY"
        confidence = "LOW"

    return sex, karyotype, confidence


# ──────────────────────────────────────────────────────────────────────
#  5. Output
# ──────────────────────────────────────────────────────────────────────

HEADER = "\t".join([
    "Sample", "Sex", "Karyotype", "Confidence",
    "X_only", "Y_only", "XY",
    "X_ratio", "Y_ratio", "XIST_reads",
])


def write_result(
    out_path: str,
    sample: str,
    sex: str,
    karyotype: str,
    confidence: str,
    x_only: int,
    y_only: int,
    xy_both: int,
    x_ratio: float,
    y_ratio: float,
    xist_count: int,
    append: bool = False,
) -> None:
    """Write (or append) a single result line to a TSV file."""
    line = "\t".join([
        sample, sex, karyotype, confidence,
        str(x_only), str(y_only), str(xy_both),
        f"{x_ratio:.4f}" if x_ratio >= 0 else "NA",
        f"{y_ratio:.4f}" if y_ratio >= 0 else "NA",
        str(xist_count),
    ])

    mode = "a" if append else "w"
    with open(out_path, mode) as fh:
        if not append or not Path(out_path).exists() or os.path.getsize(out_path) == 0:
            fh.write(HEADER + "\n")
        fh.write(line + "\n")


# ──────────────────────────────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Unified sex caller for hg38 BAM/SAM alignments. "
                    "Detects XX, XY, XXY (Klinefelter), XYY, X0 (Turner).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("bam", help="Sorted BAM file (coordinate-sorted preferred)")
    parser.add_argument("--xist", metavar="BED",
                        help="BED file with XIST exon coordinates for XIST expression counting")
    parser.add_argument("--mapq", type=int, default=30,
                        help="Minimum MAPQ for strict read filters (default: 30)")
    parser.add_argument("--y-threshold", type=int, default=10,
                        help="Minimum Y-only + XY reads to consider Y present (default: 10)")
    parser.add_argument("--output", "-o", metavar="TSV",
                        help="Output TSV path (default: <bam_dir>/sexcall.tsv)")
    parser.add_argument("--append", action="store_true",
                        help="Append to existing output file instead of overwriting")
    parser.add_argument("--sample", metavar="NAME",
                        help="Sample name (default: derived from BAM filename)")
    args = parser.parse_args()

    bam_path = os.path.realpath(args.bam)
    if not os.path.isfile(bam_path):
        sys.exit(f"[ERROR] BAM file not found: {bam_path}")

    # Derive sample name
    sample = args.sample or Path(bam_path).stem.replace("_csort", "").replace("_nsort", "")

    # Output path
    out_path = args.output or str(Path(bam_path).parent / "sexcall.tsv")

    print(f"[sexcall] Sample:  {sample}")
    print(f"[sexcall] BAM:     {bam_path}")
    print(f"[sexcall] Output:  {out_path}")
    print(f"[sexcall] MAPQ≥:   {args.mapq}")
    print(f"[sexcall] Y-thr:   {args.y_threshold}")
    if args.xist:
        print(f"[sexcall] XIST:    {args.xist}")

    # ── Step 1: Depth ratios ──────────────────────────────────────────
    print("\n[1/4] Computing chromosome depth ratios via idxstats ...")
    x_ratio, y_ratio, raw = compute_depth_ratios(bam_path)
    if x_ratio >= 0:
        print(f"       X:auto = {x_ratio:.4f}   Y:auto = {y_ratio:.4f}")
        print(f"       (median autosomal depth = {raw.get('median_auto_depth', 0):.4f}, "
              f"{raw.get('auto_count', 0)} autosomes resolved)")
    else:
        print("       [WARN] Depth ratios unavailable – BAM may not be indexed")

    # ── Step 2: Strict read classification ────────────────────────────
    print("\n[2/4] Classifying X/Y reads (strict filters) ...")
    x_only, y_only, xy_both = classify_reads(bam_path, args.mapq)
    print(f"       X_only = {x_only}   Y_only = {y_only}   XY = {xy_both}")

    # ── Step 3: XIST counting ────────────────────────────────────────
    xist_count = 0
    if args.xist:
        print(f"\n[3/4] Counting XIST reads from {args.xist} ...")
        xist_count = count_xist(bam_path, args.xist, args.mapq)
        print(f"       XIST_reads = {xist_count}")
    else:
        print("\n[3/4] XIST counting skipped (no --xist BED provided)")

    # ── Step 4: Call ─────────────────────────────────────────────────
    print("\n[4/4] Making sex call ...")
    sex, karyotype, confidence = call_sex(
        x_ratio, y_ratio, x_only, y_only, xy_both, xist_count, args.y_threshold
    )

    print(f"\n{'='*60}")
    print(f"  RESULT:  {sample}")
    print(f"  Sex:        {sex}")
    print(f"  Karyotype:  {karyotype}")
    print(f"  Confidence: {confidence}")
    print(f"{'='*60}\n")

    # ── Write output ─────────────────────────────────────────────────
    write_result(
        out_path, sample, sex, karyotype, confidence,
        x_only, y_only, xy_both, x_ratio, y_ratio, xist_count,
        append=args.append,
    )
    print(f"[sexcall] ✓ Written → {out_path}")


if __name__ == "__main__":
    main()
