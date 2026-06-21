# Analysis Notebooks

Jupyter notebooks used for the final analysis, figure generation, and report assembly.

## Notebooks

| Notebook | Description |
|----------|-------------|
| `results.ipynb` | Technical results: count-profile stability, read survival, alignment statistics, computational cost |
| `bioResults.ipynb` | Biological concordance: LOSO DE/GSEA analysis, sample-specific trimming classification, QC model evaluation |
| `methodology.ipynb` | Methodology documentation and supplementary methods |
| `bioMethodology.ipynb` | Biological analysis methodology details |

## Usage

These notebooks were developed with Python 3.14.3 and require the following key packages:
- pandas, numpy, scipy
- matplotlib, seaborn
- pydeseq2 (for differential expression)
- scikit-learn (for QC model)

The notebooks read aggregated data from the `results/` directory and per-project data from the `data/` directory.
