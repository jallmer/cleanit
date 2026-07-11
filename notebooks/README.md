# Analysis Notebooks

## Manuscript Notebooks

The `manuscript/` directory contains 11 focused notebooks, one per manuscript figure or analysis section. All are tracked on GitHub.

| Notebook | Description | Runs without Zenodo? |
|----------|-------------|---------------------|
| `00_bioproject_cohort_and_supplementary_table` | BioProject cohort overview and Supplementary Table 1 | No |
| `01_pipeline_workflow_and_tooling` | Pipeline workflow diagram and tooling summary | No |
| `02_computational_cost_and_throughput` | Computational cost and throughput analysis | No |
| `03_fastqc_quality_distribution` | FastQC quality score distributions | No |
| `04_count_profile_stability` | Count-profile stability (Pearson, Spearman, JSD) | No |
| `05_read_survival` | Read survival across trimming stringencies | No |
| `06_quality_sensitivity_to_trimming` | Quality sensitivity to trimming intensity | No |
| **`07_gene_level_concordance`** | **Gene-level biological concordance (LOSO)** | **Yes** |
| **`08_pathway_level_concordance`** | **Pathway-level concordance (GSEA + DGE)** | **Yes** |
| `09_sample_specific_trimming_classification` | Sample-specific trimming classification | No |
| `10_qc_prediction_of_trimming_need` | QC-based prediction of trimming need | No |

Notebooks **07** and **08** use only the aggregated summary tables (`results/biological/bio_concordance_binary.tsv` and `results/technical/per_srr_quality.tsv`) that are tracked on GitHub, so they can be executed immediately after cloning.

All other notebooks require the full per-project data from the [Zenodo archives](../README.md#zenodo-data).

## Dependencies

The notebooks were developed with Python 3.14.3 and mainly rely on:

- pandas, numpy, scipy
- matplotlib, seaborn
- scikit-learn
- pydeseq2 (for notebooks that regenerate DE/GSEA outputs)
