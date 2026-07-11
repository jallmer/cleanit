# Manuscript Notebook Series

This directory contains the focused notebook series for the current manuscript. Each notebook is tied to one figure, table, or larger quantitative result block.

## Notebook map

| Notebook | Manuscript target |
|----------|-------------------|
| `00_bioproject_cohort_and_supplementary_table.ipynb` | Methodology / Data Selection, Methodology / Filtering of the list of BioProjects, Supplementary Table S1 |
| `01_pipeline_workflow_and_tooling.ipynb` | Figure Workflow, Table Tools, Methodology / Data processing pipeline |
| `02_computational_cost_and_throughput.ipynb` | Methodology / Timing and Throughput Method, Results / Computational cost, Figure the computational cost |
| `03_fastqc_quality_distribution.ipynb` | Results / Quality assessment, Figure: quality distribution |
| `04_count_profile_stability.ipynb` | Methodology / Count-Profile Stability Evaluation, Results / Impact of trimming on the technical level, Figure Technical distribution (Spearman), Table count stability 1 |
| `05_read_survival.ipynb` | Methodology / Read Retention After Trimming, Results / Table read survival |
| `06_quality_sensitivity_to_trimming.ipynb` | Methodology / Integrated Quality Sensitivity Analysis, Results / Figure Spearman correlation quality plot |
| `07_gene_level_concordance.ipynb` | Results / Impact of trimming on the biological interpretation / Gene-level concordance, Figure gene level concordance |
| `08_pathway_level_concordance.ipynb` | Results / Impact of trimming on the biological interpretation / Pathway-level concordance, Figure pathway level concordance, Whole-project concordance references inside the pathway discussion |
| `09_sample_specific_trimming_classification.ipynb` | Results / Sample-specific trimming classification |
| `10_qc_prediction_of_trimming_need.ipynb` | Results / Classification of the need to trim from FastQC quality metrics, QC model figures and feature-overlap discussion |

## Design rules

- Each notebook documents its primary inputs, rebuild scripts, analysis method, and visual outputs.
- Each notebook prefers the full local `data/` working set when it exists and falls back to `share/data/` otherwise.
- The older broad notebooks in `share/notebooks/` are left in place as legacy companions until they are explicitly archived.
