# FastQC Reports

Raw FastQC quality reports for untrimmed reads.

## File Structure

Per SRR, each read direction has:
- `<SRR>_<read>_fastqc/` — extracted FastQC report directory
- `<SRR>_<read>_fastqc.html` — HTML report
- `<SRR>_<read>_fastqc.zip` — zipped report archive

Where `<read>` is `1` or `2` for paired-end data.

## Example

See `PRJNA1014965/` for sample files. Full data available on Zenodo.
