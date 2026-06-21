#!/bin/bash
# Submit one SLURM job per missing project for DESeq2 + GSEA.
#
# Usage:
#   bash scripts/submit_all_deseq2_gsea.sh

ANALYSIS="/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis"
PYTHON="/scratch/hpc-prf-omiks/ja/miniconda3/envs/omiks_pipeline/bin/python"

cd "$ANALYSIS"
mkdir -p logs

# Get list of missing projects from the Python script
MISSING=$($PYTHON -c "
from pathlib import Path
ANALYSIS = Path('$ANALYSIS')
META_DIR = ANALYSIS / 'deseq2_metadata' / 'sample_sheets'
COUNTS_DIR = ANALYSIS.parent / 'flattened_counts'
CONC_DIR = ANALYSIS / 'concordance'
for sheet in sorted(META_DIR.glob('PRJNA*.csv')):
    proj = sheet.stem
    if not (COUNTS_DIR / proj).exists():
        continue
    if (CONC_DIR / proj / 'de_U.tsv').exists() and (CONC_DIR / proj / 'gsea_U.tsv').exists():
        continue
    print(proj)
")

if [ -z "$MISSING" ]; then
    echo "All projects already have DE/GSEA results. Nothing to submit."
    exit 0
fi

COUNT=0
for PROJ in $MISSING; do
    sbatch --job-name="$PROJ" scripts/submit_deseq2_gsea.sh "$PROJ"
    COUNT=$((COUNT + 1))
done

echo ""
echo "Submitted $COUNT jobs."
