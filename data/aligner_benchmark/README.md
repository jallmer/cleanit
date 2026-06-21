# Aligner Benchmark: Cleaning × Alignment Speed

## Purpose

Independent benchmark to test whether read trimming affects alignment speed and to compare three aligners (Bowtie2, STAR, HISAT2) on the same inputs.

## Design

- **18 paired-end human RNA-seq samples** (Zika project, 0.9–3.9 GB per sample)
- **4 trimming modes**: untrimmed, adapter-only, P5, P10
- **3 aligners**: Bowtie2 (end-to-end), STAR, HISAT2
- **216 total alignment jobs**, all completed successfully
- **36 threads** per job, alignment output discarded (timing-only benchmark)

## Result

Alignment wall-clock minutes per GB of aligned FASTQ input (36 cores, n = 18):

|  | untrimmed | adapter only | P5 | P10 |
|---|---:|---:|---:|---:|
| **Bowtie2** median | 1.2 | 1.4 | 1.4 | 1.3 |
| **Bowtie2** mean | 1.5 | 1.6 | 1.6 | 1.6 |
| **STAR** median | 0.5 | 0.5 | 0.5 | 0.5 |
| **STAR** mean | 0.5 | 0.5 | 0.5 | 0.6 |
| **HISAT2** median | 1.1 | 1.1 | 1.1 | 1.1 |
| **HISAT2** mean | 2.3 | 2.3 | 2.3 | 2.3 |

**Conclusion**: Trimming does not reduce alignment time for any aligner. STAR is roughly twice as fast as Bowtie2 and HISAT2.

## Files

| File | Description |
|---|---|
| `core_hours_per_gb_per_run.tsv` | Per-SRR per-aligner per-mode results (wall time, threads, input size, alignment rate) |
| `selected_18_fastq_human.tsv` | Sample manifest (SRR IDs, conditions, FASTQ sizes, paths) |
| `run_one_fastq_aligner_mode.sh` | Slurm job script: trims + aligns one sample × one mode × one aligner |
| `summarize_aligner_benchmark.py` | Aggregation script |

## Reproducing

The `core_hours_per_gb_per_run.tsv` file contains all raw timing data. To regenerate the summary table:

```python
import pandas as pd

df = pd.read_csv("core_hours_per_gb_per_run.tsv", sep="\t")
df["wall_min_per_gb"] = (df["duration_sec"] / 60) / df["aligned_input_gb"]

for stat in ["median", "mean"]:
    pivot = df.pivot_table(
        index="mode", columns="aligner",
        values="wall_min_per_gb", aggfunc=stat
    )
    print(f"\n{stat.title()} wall-minutes per GB:")
    print(pivot.round(1).to_string())
```
