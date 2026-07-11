# Aggregated Results

Pre-computed aggregated result tables from the analysis. These stay on GitHub because they directly support the manuscript and are small enough to version here.

## `technical/`

| File | Description |
|------|-------------|
| `per_srr_eval.tsv` | Per-SRR count-profile stability metrics (Pearson, Spearman, JSD) across trimming modes |
| `per_srr_quality.tsv` | Per-SRR FastQC quality metrics (average and terminal quality) |
| `bowtie2_alignment_stats.tsv` | Bowtie2 alignment statistics per SRR and trimming mode |
| `trimmomatic_detail.tsv` | Detailed Trimmomatic read survival per SRR and trimming mode |
| `throughput_detail.tsv` | Per-stage processing times (FastQC, trimming, alignment, counting) |
| `trimming_classification.tsv` | Full per-SRR trimming classification with concordance metrics |

## `biological/`

| File | Description |
|------|-------------|
| `bio_concordance_binary.tsv` | LOSO biological concordance metrics (DGE + GSEA: Spearman ρ, Jaccard, direction concordance) per sample × trimming mode, using true binary control/treatment groupings |
| `classification_summary.tsv` | Summary of sample-specific trimming classification (helpful/neutral/harmful counts per method) |
| `qc_model_predictions.tsv` | Random Forest QC model predictions and feature importances |
| `loso_binary/` | Full LOSO analysis outputs |
| `loso_binary/deseq2_sample_annotation.tsv` | Binary group annotations for all 41 BioProjects |
| `loso_binary/results/PRJNA1014965_concordance.tsv` | Example per-project concordance result (full set on [Zenodo](https://doi.org/10.5281/zenodo.21296496)) |
| `loso_binary/scripts/06b_run_loso_binary.py` | Script to reproduce the LOSO analysis |

## `supplementary/`

| File | Description |
|------|-------------|
| `Supplementary_Table_1.csv` | Supplementary Table 1 for the manuscript |
