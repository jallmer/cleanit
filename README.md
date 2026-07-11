# CleanIt - Shared Materials

**On the need to trim RNA-seq data for transcriptomic data analysis**

Fabian Bräuer and Jens Allmer

---

## Overview

This directory is the publication-facing `share` repo for CleanIt. GitHub keeps the manuscript-facing subset here: aggregated results, reference resources, and a few example flattened files. Full scripts, notebooks, SQLite state, and the complete flattened HPC outputs are distributed via Zenodo and may exist locally in this directory without being synced to GitHub.

## Distribution

| Channel | Contents |
|---------|----------|
| **GitHub** (this repository) | READMEs, aggregated result tables, reference resources, and a few curated example files from `PRJNA1014965` |
| **Zenodo** ([DOI: 10.5281/zenodo.20787459](https://doi.org/10.5281/zenodo.20787459)) | Full scripts, notebooks, SQLite state, and the complete flattened HPC outputs |

## Directory Layout

```
share/
├── results/               Aggregated analysis outputs tracked on GitHub
│   ├── technical/         Count-profile stability, alignment stats, trimming stats
│   ├── biological/        LOSO concordance, trimming classification, QC model
│   └── supplementary/     Supplementary tables for the manuscript
├── resources/             Reference inputs tracked on GitHub
├── data/                  Example flattened files on GitHub; full bundles on Zenodo
│   ├── flattened_counts/
│   ├── flattened_bowtie2_stats/
│   ├── flattened_fastqc_raw/
│   ├── flattened_trimmomatic_stats/
│   ├── flattened_timings/
│   └── concordance/
├── scripts/               Local/Zenodo-only operational code
├── notebooks/             Local/Zenodo-only notebooks
└── srr_queue.db           Local/Zenodo-only SQLite pipeline database
```

## Reproducing the Analysis

1. Download the Zenodo release and unpack it into `share/` to populate the full `data/`, `scripts/`, `notebooks/`, and `srr_queue.db` payloads.
2. Use the GitHub-tracked `results/` and `resources/` trees for manuscript-facing tables and reference inputs.
3. See [`data/README.md`](./data/README.md), [`scripts/README.md`](./scripts/README.md), and [`notebooks/README.md`](./notebooks/README.md) for the expected local layout after unpacking Zenodo.

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
