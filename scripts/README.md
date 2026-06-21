# Scripts

This directory contains all scripts used in the CleanIt analysis, organized by function.

## Subdirectories

### `pipeline/`
The HPC processing pipeline scripts. The main entry point is `v2_master_node.sh`, which orchestrates SRA conversion, FastQC, Trimmomatic trimming (six modes), Bowtie2 alignment, and featureCounts counting for each SRR.

### `data_acquisition/`
Scripts for downloading SRR data from NCBI and EBI endpoints (SRA, AWS, GCP, ENA HTTP, ENA Aspera, DDBJ). The pipeline uses endpoint fallback to handle download failures.

### `analysis/`
Numbered analysis scripts that process HPC outputs into the final results. Run in order from `00_setup_env.sh` through the numbered scripts. Key scripts:
- `02_evaluate_count_matrices.py` — count-profile stability (Pearson, Spearman, JSD)
- `06_run_de_gsea.py` — LOSO differential expression and GSEA
- `07_classify_trimming.py` — sample-specific trimming classification
- `08_fit_qc_model.py` — Random Forest QC prediction model
- `11_compute_throughput.py` — computational cost analysis

### `utilities/`
Helper scripts for ancillary analyses (sex calling, viral coverage profiling, archive management).
