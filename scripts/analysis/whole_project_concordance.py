#!/usr/bin/env python
"""
Compute whole-project biological concordance metrics.

For each eligible BioProject that has DESeq2 and GSEA results for all trimming
modes, compare trimmed results against the untrimmed reference at full sample
power. This fills the analytical gap between the technical count-profile
analysis and the LOSO robustness evaluation.

Reads:
  concordance/<PROJECT>/de_<M>.tsv   — DESeq2 results per method
  concordance/<PROJECT>/gsea_<M>.tsv — GSEA prerank results per method

Produces:
  concordance/whole_project_concordance.tsv — per-project, per-method metrics

Usage:
    python scripts/whole_project_concordance.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr

# ── Paths ──
ANALYSIS = Path(__file__).resolve().parent.parent
CONCORDANCE_DIR = ANALYSIS / "concordance"
METHODS = ["A", "P5", "P10", "P20", "P35"]

DEG_PADJ = 0.05
DEG_LFC  = 1.0
GSEA_FDR = 0.25


def load_de(proj_dir, method):
    f = proj_dir / f"de_{method}.tsv"
    return pd.read_csv(f, sep="\t", index_col=0) if f.exists() else None


def load_gsea(proj_dir, method):
    f = proj_dir / f"gsea_{method}.tsv"
    return pd.read_csv(f, sep="\t", index_col=0) if f.exists() else None


def get_degs(de_df):
    if de_df is None:
        return set()
    mask = (de_df["padj"].fillna(1) < DEG_PADJ) & (de_df["log2FoldChange"].abs() > DEG_LFC)
    return set(de_df.index[mask])


def get_sig_pathways(gsea_df):
    if gsea_df is None:
        return set()
    return set(gsea_df.loc[gsea_df["FDR q-val"].fillna(1) < GSEA_FDR, "Term"])


def jaccard(a, b):
    return len(a & b) / len(a | b) if len(a | b) > 0 else np.nan


def discover_projects():
    projects = []
    for d in sorted(CONCORDANCE_DIR.iterdir()):
        if d.is_dir() and (d / "de_U.tsv").exists() and (d / "gsea_U.tsv").exists():
            projects.append(d.name)
    return projects


def main():
    projects = discover_projects()
    print(f"Found {len(projects)} projects with whole-project DE + GSEA results")

    rows = []
    for proj in projects:
        proj_dir = CONCORDANCE_DIR / proj
        de_u   = load_de(proj_dir, "U")
        gsea_u = load_gsea(proj_dir, "U")
        degs_u = get_degs(de_u)
        paths_u = get_sig_pathways(gsea_u)

        for m in METHODS:
            de_m   = load_de(proj_dir, m)
            gsea_m = load_gsea(proj_dir, m)
            degs_m = get_degs(de_m)
            paths_m = get_sig_pathways(gsea_m)

            # DEG metrics
            rho_lfc = np.nan
            if de_u is not None and de_m is not None:
                common = de_u.index.intersection(de_m.index)
                lfc_u = de_u.loc[common, "log2FoldChange"].dropna()
                lfc_m = de_m.loc[common, "log2FoldChange"].dropna()
                shared = lfc_u.index.intersection(lfc_m.index)
                if len(shared) > 10:
                    rho_lfc, _ = spearmanr(lfc_u[shared], lfc_m[shared])

            # Top-50 DEG rank stability
            rho_rank = np.nan
            if de_u is not None and de_m is not None:
                stat_u = de_u["stat"].dropna().abs().sort_values(ascending=False)
                stat_m = de_m["stat"].dropna().abs()
                top50 = stat_u.head(50).index.intersection(stat_m.index)
                if len(top50) > 5:
                    rho_rank, _ = spearmanr(stat_u[top50], stat_m[top50])

            # GSEA metrics
            rho_nes = np.nan
            dir_conc = np.nan
            if gsea_u is not None and gsea_m is not None:
                merged = gsea_u[["Term", "NES"]].merge(
                    gsea_m[["Term", "NES"]], on="Term", suffixes=("_U", "_M"))
                merged["NES_U"] = pd.to_numeric(merged["NES_U"], errors="coerce")
                merged["NES_M"] = pd.to_numeric(merged["NES_M"], errors="coerce")
                merged = merged.dropna(subset=["NES_U", "NES_M"])
                if len(merged) > 3:
                    rho_nes, _ = spearmanr(merged["NES_U"], merged["NES_M"])
                top20 = merged.assign(
                    absNES=merged["NES_U"].abs()
                ).nlargest(20, "absNES")
                if len(top20) > 0:
                    same = (np.sign(top20["NES_U"]) == np.sign(top20["NES_M"])).sum()
                    dir_conc = same / len(top20)

            rows.append({
                "project": proj, "method": m,
                "n_deg_untrimmed": len(degs_u), "n_deg_trimmed": len(degs_m),
                "deg_jaccard": jaccard(degs_u, degs_m),
                "deg_gained": len(degs_m - degs_u),
                "deg_lost": len(degs_u - degs_m),
                "rho_log2fc": rho_lfc,
                "rho_top50_rank": rho_rank,
                "n_path_untrimmed": len(paths_u), "n_path_trimmed": len(paths_m),
                "path_jaccard": jaccard(paths_u, paths_m),
                "rho_nes": rho_nes,
                "dir_concordance": dir_conc,
            })

    summary = pd.DataFrame(rows)
    out = CONCORDANCE_DIR / "whole_project_concordance.tsv"
    summary.to_csv(out, sep="\t", index=False)
    print(f"Saved → {out}")
    print(f"  {len(summary)} rows ({len(projects)} projects × {len(METHODS)} methods)")


if __name__ == "__main__":
    main()
