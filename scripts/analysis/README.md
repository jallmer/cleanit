# Analysis Scripts

Post-processing analysis scripts. The numbered scripts form the main analysis chain; additional scripts handle specific sub-analyses.

## Main Analysis Chain

| Script | Purpose |
|--------|---------|
| `00_setup_env.sh` | Environment setup |
| `01_extract_fastqc_stats.py` / `.sh` | Extract quality metrics from FastQC outputs |
| `02_evaluate_count_matrices.py` | Count-profile stability (Pearson, Spearman, cosine, JSD) |
| `03_extract_trimmomatic_stats.sh` | Parse Trimmomatic logs for read survival |
| `04_summarize_results.py` | Aggregate per-SRR evaluation metrics |
| `05_run_all.sh` | Master runner for the analysis chain |
| `06_run_de_gsea.py` | LOSO differential expression and GSEA concordance |
| `07_classify_trimming.py` | Sample-specific trimming classification |
| `08_fit_qc_model.py` | Train Random Forest QC prediction model |
| `09_generate_reports.py` | Generate result reports and figures |
| `10_fetch_metadata.py` | Fetch BioProject metadata from NCBI |
| `11_compute_throughput.py` | Computational cost analysis |

## Extended Analysis Scripts

| Script | Purpose |
|--------|---------|
| `12_backfill_flattened_stats_from_logs.py` | Recover stats from pipeline logs |
| `13_build_sample_feature_table.py` | Build feature table for QC model |
| `13_merge_fb_results.py` | Merge featureBranch results |
| `14_merge_concordance.py` | Merge LOSO concordance tables |
| `14_merge_fb_bowtie_results.py` | Merge featureBranch Bowtie2 stats |
| `14_regenerate_diagnostics.py` | Regenerate diagnostic plots |
| `15_deduplicate_alignment_stats.py` | Deduplicate alignment statistics |
| `15_plot_quality_balanced_correlation.py` | Quality-balanced correlation plots |
| `16_integrate_fb_counts.py` | Integrate featureBranch counts |
| `16_plot_tail_quality_by_read_length.py` | Tail-quality vs read-length plots |
| `17_build_deseq2_sample_sheets.py` | Build DESeq2 sample sheets from metadata |
| `17_extract_fb_counting_timing.py` | Extract counting timings |
| `18_integrate_fb_timings.py` | Integrate timing data |

## LOSO / DE-GSEA Submission

| Script | Purpose |
|--------|---------|
| `loso_job.sh` | Single LOSO job for Slurm |
| `submit_all_loso.sh` | Submit all LOSO jobs |
| `submit_dynamic_loso.py` | Dynamic LOSO submission |
| `run_bio_concordance.py` | Compute biological concordance |
| `run_qc_fitting.sh` | End-to-end QC model fitting chain |

## Report / Notebook Helpers

| Script | Purpose |
|--------|---------|
| `add_final_sections.py` | Add sections to results notebook |
| `append_confusion_matrix.py` | Add confusion matrix to report |
| `append_notebook.py` | Append cells to notebook |
| `build_package.py` | Package results |
| `generate_corrected_table.py` | Generate corrected result tables |

## Whole-Project and Classification Scripts

| Script | Purpose |
|--------|---------|
| `aggregate_loso_concordance.py` | Aggregate LOSO concordance results |
| `generate_classification.py` | Generate trimming classifications |
| `whole_project_concordance.py` | Whole-project (non-LOSO) concordance |
| `run_deseq2_gsea.py` | Run DESeq2 + GSEA for a single project |
| `submit_all_deseq2_gsea.sh` / `submit_deseq2_gsea.sh` | Slurm submission wrappers |
| `retrain_qc_model.py` | Retrain QC model with updated data |
| `results_technical.py` | Technical results table generation |
| `rebuild_quality_table.py` | Rebuild quality summary table |

## Notebook Update Helpers

| Script | Purpose |
|--------|---------|
| `add_11abc.py` / `add_9c_9d.py` | Add specific notebook sections |
| `add_figures_to_notebook.py` | Insert figures into notebooks |
| `update_cell22.py` / `update_notebook_9ab.py` / `update_notebook_tstar.py` | Update specific notebook cells |

## Technical Archive Scripts

| Script | Purpose |
|--------|---------|
| `archive_and_prune_filtered_srrs.py` | Archive filtered SRR outputs |
| `archive_pruned_final_fastqc_files.py` | Archive pruned FastQC files |
| `build_final_fastqc_bundle.py` | Build final FastQC data bundle |
| `build_final_trimming_bundle.py` | Build final trimming data bundle |
| `extract_final_fastqc_srr_metrics.py` | Extract per-SRR FastQC metrics |
| `prune_filtered_srr_rows_from_tsvs.py` | Remove filtered SRRs from tables |
| `prune_final_fastqc_to_db_srrs.py` | Prune FastQC to database SRRs |
