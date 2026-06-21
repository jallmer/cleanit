# Data Acquisition Scripts

Scripts for downloading sequencing data from multiple endpoints. The pipeline falls back across endpoints when downloads fail.

## Endpoints

| Script | Source | Protocol |
|--------|--------|----------|
| `fetch_ncbi.sh` | NCBI SRA | HTTPS |
| `fetch_ncbi_cloud.sh` | NCBI via cloud | AWS/GCP URLs |
| `fetch_aws.sh` | AWS S3 | HTTPS |
| `fetch_gcp.sh` | Google Cloud | HTTPS |
| `fetch_ena.sh` | ENA | Generic |
| `fetch_ena_http.sh` | ENA | HTTP |
| `fetch_ena_aspera.sh` | ENA | Aspera |
| `fetch_ddbj.sh` | DDBJ | HTTPS |
| `fetch_ddbj_sra.sh` | DDBJ | SRA format |

## Dispatchers

| Script | Role |
|--------|------|
| `v2_data_fetcher.sh` | Main data fetcher with endpoint fallback |
| `v2_fetch_dispatcher.sh` | Dispatches fetch jobs across endpoints |
| `v2_fetch_worker.sh` | Worker process for parallel downloads |
| `v2_slurm_fetch.sh` | Slurm-submitted fetch job |
| `fetch_data.sh` | Earlier fetch wrapper |
