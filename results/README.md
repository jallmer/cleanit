# Aggregated Results

Pre-computed aggregated result tables from the analysis. These are small enough to include directly in the GitHub repository.

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
| `bio_concordance.tsv` | LOSO biological concordance metrics (Spearman ρ, Jaccard, direction concordance) per sample × trimming mode |
| `whole_project_concordance.tsv` | Whole-project (non-LOSO) concordance: DEG counts, Jaccard overlap, effect-size correlation |
| `classification_summary.tsv` | Summary of sample-specific trimming classification (helpful/neutral/harmful counts per method) |
| `concordance_summary_by_method.tsv` | Aggregated concordance statistics per trimming method |
| `qc_model_predictions.tsv` | Random Forest QC model predictions and feature importances |

## `supplementary/`

| File | Description |
|------|-------------|
| `Supplementary_Table_1.csv` | Supplementary Table 1 for the manuscript |
