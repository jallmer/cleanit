#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────
# extract_xist_bed.sh – Generates an XIST BED file from a GTF/GFF
# ─────────────────────────────────────────────────────────────────────
# This script extracts all exon coordinates for the XIST gene from an
# input GTF file and converts them into a 0-indexed BED format, which
# is required by SexCall for XXY vs XY disambiguation.
# ─────────────────────────────────────────────────────────────────────

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <fasta_file> <gtf_file> [output_bed]"
    echo "Example: $0 GRCh38.fasta gencode.v44.annotation.gtf XIST_exons.GRCh38.bed"
    exit 1
fi

FASTA_FILE="$1"
GTF_FILE="$2"

if [[ $# -ge 3 ]]; then
    OUTPUT="$3"
else
    BASENAME=$(basename "${FASTA_FILE}" | sed -E 's/\.(fna|fa|fasta)(\.gz)?$//i')
    OUTPUT="XIST_exons.${BASENAME}.bed"
fi

if [[ ! -f "$GTF_FILE" ]]; then
    echo "[ERROR] GTF file not found: $GTF_FILE" >&2
    exit 1
fi

echo "🔹 Extracting XIST exons from $GTF_FILE ..."

# Set up command to read either plain text or gzipped GTF
CMD="cat"
if [[ "$GTF_FILE" == *.gz ]]; then
    CMD="zcat"
fi

# We look for lines where feature (col 3) is "exon" and the attributes (col 9) 
# contain either gene_name "XIST" (GENCODE/Ensembl) or gene "XIST" (RefSeq).
# GTF coordinates are 1-based, inclusive. BED coordinates are 0-based, exclusive.
# Therefore, BED start = GTF start - 1.

$CMD "$GTF_FILE" | awk -F'\t' '
BEGIN { OFS="\t" }
# Skip header lines
/^#/ { next }
# Match feature == "exon" and XIST in attributes
$3 == "exon" && ($9 ~ /gene_name "XIST"/ || $9 ~ /gene "XIST"/) {
    start = $4 - 1
    end = $5
    # Output: chrom, start, end, name, score, strand
    print $1, start, end, "XIST", ".", $7
}
' > "$OUTPUT"

COUNT=$(wc -l < "$OUTPUT" | awk '{print $1}')

if [[ "$COUNT" -gt 0 ]]; then
    echo "✅ Successfully wrote $COUNT XIST exons to $OUTPUT"
else
    echo "⚠️  [WARNING] No XIST exons found in the GTF file! Ensure this is a human hg38/GRCh38 annotation."
    # Optionally, we leave the empty file or remove it. We'll leave it in case it's intentionally empty for a non-human run.
fi
