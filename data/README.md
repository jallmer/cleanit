# Per-Sample HPC Data

This directory contains **sample reference files** from one example BioProject (`PRJNA1014965`) so that users can see the file naming conventions and directory layout.

**The full dataset is available on Zenodo: [DOI: 10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)**

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
├── flattened_timings/          Per-stage pipeline timings
│   └── <BioProject>/
│       └── <SRR>_timings.tsv
└── concordance/                LOSO DE/GSEA concordance results
    └── <BioProject>/
        └── (per-sample concordance files)
```

## Trimming Modes

The `<mode>` in count file names refers to:
- `adapter` — adapter-only trimming (ILLUMINACLIP only)
- `P5` — adapter + PHRED 5 quality trimming
- `P10` — adapter + PHRED 10 quality trimming
- `P20` — adapter + PHRED 20 quality trimming
- `P35` — adapter + PHRED 35 quality trimming

Untrimmed counts use the raw reads directly (no Trimmomatic processing).

## Total Data Volume

The full dataset spans 48 BioProjects with ~1,155 SRRs, each processed under 6 trimming conditions.
