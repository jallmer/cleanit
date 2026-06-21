# Pipeline Scripts

HPC processing pipeline for the CleanIt project. These scripts run on a Slurm-managed cluster.

## Main Entry Point

**`v2_master_node.sh`** — Orchestrates the full per-SRR pipeline:
1. SRA data download and FASTQ conversion
2. FastQC quality assessment
3. Trimmomatic trimming (6 modes: untrimmed, adapter-only, P5, P10, P20, P35)
4. Bowtie2 alignment (end-to-end mode, no soft clipping)
5. Samtools BAM conversion and sorting
6. featureCounts gene-level counting

## Supporting Scripts

| Script | Role |
|--------|------|
| `v2_trim.sh` | Trimmomatic wrapper for all trimming modes |
| `v2_align.sh` | Bowtie2 alignment wrapper |
| `v2_count.sh` | featureCounts wrapper |
| `v2_supervisor.sh` | Job monitoring and resubmission |
| `v2_submitter.sh` / `v2_pure_submitter.sh` | Slurm job submission helpers |
| `v2_submit_supervisor.sh` | Submits the supervisor job |
| `v2_macro_runner.sh` | Batch runner for multiple projects |
| `pipeline.sh` | Earlier single-script pipeline version |
| `reconcile_db.py` | Database reconciliation utility |
| `resubmit_pipelines.sh` | Resubmit failed pipeline jobs |
| `start_missing_pipeline.sh` | Identify and start unprocessed SRRs |
