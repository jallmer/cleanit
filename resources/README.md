# Resources

Reference files, adapter sequences, project lists, and sample metadata used throughout the analysis. These remain on GitHub because they are small manuscript-facing inputs.

## Files

| File | Description |
|------|-------------|
| `adapter_sequences.fasta` | Illumina adapter sequences used by Trimmomatic (from the [Trimmomatic GitHub repository](https://github.com/usadellab/Trimmomatic)) |
| `project_IDs_filtered.tsv` | 55 BioProject accessions passing the initial size filter (10–50 SRRs) |
| `project_IDs_biological_eval.tsv` | 42 BioProjects used for the biological (DE/GSEA) evaluation |
| `MSigDB_Hallmark_2020.gmt` | MSigDB Hallmark 2020 gene sets used for GSEA |
| `Supplementary_Table_BioProjects.csv` | Supplementary table with BioProject metadata |

## Subdirectory

### `project_metadata/`
Per-project CSV files mapping SRR accessions to experimental groups (conditions/classes) for differential expression analysis. Each file has one row per sample with columns for the SRR ID and group label.
