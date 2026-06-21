# featureCounts Count Files

Gene-level count matrices produced by featureCounts (Subread v2.1.1) against the GRCh38.p14 annotation.

## File Naming

- `<SRR>_trimmomatic_<mode>_fC.txt.gz` — count table (gene ID, count)
- `<SRR>_trimmomatic_<mode>_fC.txt.summary.gz` — featureCounts assignment summary

## Modes

`adapter`, `P5`, `P10`, `P20`, `P35`. Untrimmed counts are separate (no `_trimmomatic_` prefix).

## Example

See `PRJNA1014965/` for sample files. Full data available on Zenodo.
