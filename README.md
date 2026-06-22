# CleanIt — Shared Materials

**On the need to trim RNA-seq data for transcriptomic data analysis**

Fabian Bräuer and Jens Allmer

---

## Overview

This directory contains all scripts, analysis notebooks, supporting resources, and aggregated results from the CleanIt project. The study evaluates whether read trimming is necessary for standard Illumina RNA-seq differential expression workflows by processing approximately one thousand samples across six trimming conditions.

## Distribution

| Channel | Contents |
|---------|----------|
| **GitHub** (this repository) | Scripts, notebooks, resources, aggregated result tables, sample data files |
| **Zenodo** ([DOI: 10.5281/zenodo.20787459](https://doi.org/10.5281/zenodo.20787459)) | Full per-sample HPC outputs: featureCounts count files, Bowtie2 alignment stats, FastQC reports, Trimmomatic logs, pipeline timings |

## Directory Layout

```
share/
├── scripts/               Pipeline and analysis scripts
│   ├── pipeline/          HPC processing pipeline (trimming, alignment, counting)
│   ├── data_acquisition/  SRR download and fetch scripts
│   ├── analysis/          Post-processing analysis (DE, GSEA, QC model, reports)
│   └── utilities/         Helper scripts (sex calling, viral profiling, archiving)
├── notebooks/             Jupyter notebooks for results and methodology
├── resources/             Adapter sequences, project lists, sample metadata, gene sets
├── results/               Aggregated analysis outputs (TSV tables)
│   ├── technical/         Count-profile stability, alignment stats, trimming stats
│   ├── biological/        LOSO concordance, trimming classification, QC model
│   └── supplementary/     Supplementary tables for the manuscript
├── data/                  Per-sample HPC outputs (sample files; full data on Zenodo)
│   ├── flattened_counts/
│   ├── flattened_bowtie2_stats/
│   ├── flattened_fastqc_raw/
│   ├── flattened_trimmomatic_stats/
│   ├── flattened_timings/
│   └── concordance/
└── srr_queue.db           SQLite database used for pipeline orchestration
```

## Reproducing the Analysis

1. Download the full per-sample data from Zenodo and place it under `data/` following the directory layout shown by the sample files.
2. See `scripts/analysis/` and the numbered scripts (`00_setup_env.sh` through `11_compute_throughput.py`) for the analysis pipeline.
3. Open the notebooks in `notebooks/` to reproduce figures and tables.

## Software Versions

| Tool | Version |
|------|---------|
| SRA Toolkit | 3.4.1 |
| FastQC | 0.12.1 |
| Trimmomatic | 0.40 |
| Bowtie2 | 2.5.4 |
| Samtools | 1.23 |
| featureCounts / Subread | 2.1.1 |
| Python | 3.14.3 |

## Citation

If you use these materials, please cite:

> Bräuer, F. and Allmer, J. (2026). On the need to trim RNA-seq data for transcriptomic data analysis. *[Journal TBD]*.
