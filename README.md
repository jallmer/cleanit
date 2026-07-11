# CleanIt - Shared Materials

**On the need to trim RNA-seq data for transcriptomic data analysis**

Fabian Bräuer and Jens Allmer

---

## Overview

This directory is the publication-facing `share` repo for CleanIt. GitHub tracks the manuscript notebooks, aggregated result tables, reference resources, and a few example flattened files from one BioProject (`PRJNA1014965`). Full per-project raw data, pipeline scripts, and bulk HPC outputs are distributed via Zenodo.

## Zenodo Data

The raw data and bulk exports are hosted across four Zenodo records. To fully reproduce all analyses locally, download and extract the corresponding datasets:

| Zenodo Record | DOI | Contents | Extract Into |
|---------------|-----|----------|-------------|
| Alignment counts | [10.5281/zenodo.21296565](https://doi.org/10.5281/zenodo.21296565) | `cleanit_counts.tar.gz` — gene-level featureCounts matrices for all BioProjects | `data/flattened_counts/` |
| LOSO analysis results | [10.5281/zenodo.21296496](https://doi.org/10.5281/zenodo.21296496) | `loso_binary_full_results.tar.gz` — per-project LOSO concordance results | `results/biological/loso_binary/` |
| Technical stats | [10.5281/zenodo.21306331](https://doi.org/10.5281/zenodo.21306331) | `cleanit_timings.tar.gz`, `cleanit_trimmomatic_stats.tar.gz`, `cleanit_bowtie2_stats.tar.gz`, `cleanit_fastqc_raw.tar.gz` | `data/` |
| Concordance tables | [10.5281/zenodo.21306381](https://doi.org/10.5281/zenodo.21306381) | `cleanit_concordance.tar.gz` — pre-computed concordance tables | `results/biological/` |

> **Note:** The manuscript notebooks (07 and 08) can run directly from the GitHub-tracked summary tables without downloading any Zenodo data. The remaining notebooks require the Zenodo data for full execution.

## Directory Layout

```
share/
├── notebooks/
│   └── manuscript/            All 11 manuscript notebooks (tracked on GitHub)
│       ├── 00_bioproject_cohort_and_supplementary_table.ipynb
│       ├── ...
│       ├── 07_gene_level_concordance.ipynb
│       ├── 08_pathway_level_concordance.ipynb
│       └── 10_qc_prediction_of_trimming_need.ipynb
├── results/                   Aggregated analysis outputs (tracked on GitHub)
│   ├── technical/             Count-profile stability, alignment stats, trimming stats
│   ├── biological/            LOSO concordance, classification, QC model predictions
│   │   ├── bio_concordance_binary.tsv
│   │   └── loso_binary/       Example result + reproduction script
│   └── supplementary/         Supplementary tables for the manuscript
├── resources/                 Reference inputs (tracked on GitHub)
├── data/                      Example flattened files on GitHub; full bundles on Zenodo
│   ├── flattened_counts/
│   ├── flattened_bowtie2_stats/
│   ├── flattened_fastqc_raw/
│   ├── flattened_trimmomatic_stats/
│   └── flattened_timings/
└── scripts/                   Pipeline and analysis scripts (Zenodo/local only)
```

## Reproducing the Analysis

1. Clone this repository.
2. Download the Zenodo archives listed above and extract each into the indicated directory.
3. The manuscript notebooks in `notebooks/manuscript/` are numbered in reading order.
   - Notebooks **07** and **08** (biological concordance) work out of the box with just the GitHub data.
   - All other notebooks require the Zenodo data to be present locally.
4. See [`results/README.md`](./results/README.md), [`data/README.md`](./data/README.md), and [`notebooks/README.md`](./notebooks/README.md) for details.

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
