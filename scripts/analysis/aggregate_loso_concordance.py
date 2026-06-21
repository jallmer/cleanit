#!/usr/bin/env python
"""
Aggregate per-project LOSO concordance TSVs into a single long-format table.

Reads every PRJNA*_concordance.tsv from the concordance/ directory, drops
failed LOSO folds (rows with NA), and melts the wide-format per-method columns
into a tidy long-format table suitable for bioResults.ipynb.

Output: concordance/bio_concordance.tsv
        Columns: project, SRR_ID, method, rho_gene, jaccard_deg,
                 rho_pathway, jaccard_pathway, dir_concordance

Usage:
    python scripts/aggregate_loso_concordance.py
"""

import pandas as pd
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parent.parent
CONCORDANCE_DIR = ANALYSIS / "concordance"

METHODS = ["U", "A", "P5", "P10", "P20", "P35"]
METRICS_MAP = {
    "rho_gene":       "rho_gene",
    "jaccard_deg":    "jaccard_deg",
    "rho_path":       "rho_pathway",
    "jaccard_path":   "jaccard_pathway",
    "dir_concordance":"dir_concordance",
}


def main():
    files = sorted(CONCORDANCE_DIR.glob("PRJNA*_concordance.tsv"))
    if not files:
        print("No PRJNA*_concordance.tsv files found.")
        return

    rows = []
    for f in files:
        if f.name in ("all_concordance.tsv", "bio_concordance.tsv"):
            continue
        try:
            df_wide = pd.read_csv(f, sep="\t")
        except Exception as e:
            print(f"  SKIP {f.name}: {e}")
            continue

        df_wide = df_wide.dropna()  # Remove failed LOSO folds

        for _, row in df_wide.iterrows():
            for m in METHODS:
                new_row = {
                    "project": row["project_id"],
                    "SRR_ID":  row["SRR_ID"],
                    "method":  m,
                }
                for wide_key, long_key in METRICS_MAP.items():
                    new_row[long_key] = row.get(f"{m}_{wide_key}", None)
                rows.append(new_row)

    df_long = pd.DataFrame(rows).dropna()

    out = CONCORDANCE_DIR / "bio_concordance.tsv"
    df_long.to_csv(out, sep="\t", index=False)

    print(f"Aggregated {len(files)} project files")
    print(f"  {df_long['SRR_ID'].nunique()} samples, "
          f"{df_long['project'].nunique()} projects")
    print(f"  {len(df_long)} total rows (5 methods × valid samples)")
    print(f"Saved → {out}")


if __name__ == "__main__":
    main()
