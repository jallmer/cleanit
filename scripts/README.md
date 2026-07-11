# Scripts

The operational scripts are distributed via Zenodo and are not part of the GitHub-tracked subset of `share/`. Keep local copies here if you are rebuilding packages or rerunning the workflow; GitHub intentionally keeps only this README.

## Subdirectories

### `pipeline/`
The end-to-end HPC pipeline. This bucket now includes SRR fetching, queue orchestration, Slurm submission, SRA conversion, FastQC, Trimmomatic trimming, Bowtie2 alignment, and featureCounts counting. Files are ordered by stage, with `00_...` reserved for the full per-SRR pipeline entry point.

### `analysis/`
Post-flattening manuscript analysis scripts. The active workflow now reads in order from `00_install_analysis_environment.sh` through `11_compute_pipeline_throughput.py`, with additional numbered helpers for LOSO aggregation, whole-project DE/GSEA, QC model validation, and sex-calling support.

Scripts that were duplicates, legacy one-offs, bundle/pruning utilities, or the separate aligner benchmark were moved out of the active path into `archive/legacy_code/share_scripts_analysis/`.

### `local/`
Repo-facing notebook and manuscript helpers. These scripts edit notebooks, assemble report sections, and generate the share-facing statistical notebook and plots from repository-relative paths.
