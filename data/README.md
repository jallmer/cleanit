# Per-Sample HPC Data

This directory keeps a small set of **sample reference files** from one example BioProject (`PRJNA1014965`) so users can see the naming conventions and directory layout on GitHub.

**The full datasets are available on Zenodo:**

- **Alignment counts:** [10.5281/zenodo.21296565](https://doi.org/10.5281/zenodo.21296565) — extract `cleanit_counts.tar.gz` into `data/flattened_counts/`
- **Technical stats:** [10.5281/zenodo.21306331](https://doi.org/10.5281/zenodo.21306331) — extract `cleanit_fastqc_raw.tar.gz`, `cleanit_bowtie2_stats.tar.gz`, `cleanit_trimmomatic_stats.tar.gz`, and `cleanit_timings.tar.gz` into the corresponding `data/flattened_*` subdirectories

After downloading from Zenodo, extract the archives into the corresponding subdirectories here. Each subdirectory is organized by BioProject accession, with per-SRR files inside.

## Directory Layout

```
data/
├── flattened_counts/           Gene-level featureCounts output
│   └── <BioProject>/
│       ├── <SRR>_trimmomatic_<mode>_fC.txt.gz
│       └── <SRR>_trimmomatic_<mode>_fC.txt.summary.gz
├── flattened_bowtie2_stats/    Bowtie2 alignment statistics
│   └── <BioProject>/
│       └── <SRR>_bowtie2_stats.tsv
├── flattened_fastqc_raw/       Raw FastQC reports
│   └── <BioProject>/
│       ├── <SRR>_<read>_fastqc/
│       ├── <SRR>_<read>_fastqc.html
│       └── <SRR>_<read>_fastqc.zip
├── flattened_trimmomatic_stats/ Trimmomatic read survival logs
│   └── <BioProject>/
│       └── <SRR>_trimmomatic_stats.tsv
└── flattened_timings/          Per-stage pipeline timings
    └── <BioProject>/
        └── <SRR>_timings.tsv
```

## Trimming Modes

The `<mode>` in count file names refers to:
- `adapter` - adapter-only trimming (ILLUMINACLIP only)
- `P5` - adapter + PHRED 5 quality trimming
- `P10` - adapter + PHRED 10 quality trimming
- `P20` - adapter + PHRED 20 quality trimming
- `P35` - adapter + PHRED 35 quality trimming

Untrimmed counts use the raw reads directly (no Trimmomatic processing).

## Scope

The GitHub copy is intentionally minimal. Use it to inspect the layout and naming scheme; use Zenodo for the full per-project and per-SRR payloads.
