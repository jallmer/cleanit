#!/usr/bin/env python
"""
Run whole-project DESeq2 + GSEA for specified BioProjects.

For each project, builds a count matrix from featureCounts output for each
trimming mode (U, A, P5, P10, P20, P35), runs PyDESeq2 differential expression,
then GSEA prerank. Results are saved to concordance/<PROJECT>/de_<M>.tsv and
gsea_<M>.tsv.

Inputs:
  - flattened_counts/<PROJECT>/  — featureCounts .gz files
  - deseq2_metadata/sample_sheets/<PROJECT>.csv  — sample sheet with Run,condition
  - MSigDB_Hallmark_2020.gmt  — gene set library

Outputs:
  - concordance/<PROJECT>/de_<M>.tsv   — DESeq2 results per method
  - concordance/<PROJECT>/gsea_<M>.tsv — GSEA prerank results per method

Usage:
    # Run for all projects missing DE/GSEA results:
    python scripts/run_deseq2_gsea.py

    # Run for specific projects:
    python scripts/run_deseq2_gsea.py PRJNA1105191 PRJNA381757
"""

import sys
import gzip
import warnings
import pandas as pd
import numpy as np
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ──
ANALYSIS    = Path(__file__).resolve().parent.parent
COUNTS_DIR  = ANALYSIS.parent / "flattened_counts"
META_DIR    = ANALYSIS / "deseq2_metadata" / "sample_sheets"
CONC_DIR    = ANALYSIS / "concordance"
GMT_FILE    = ANALYSIS / "MSigDB_Hallmark_2020.gmt"

METHOD_MAP = {
    "U":   "untrmd_{srr}_fC.txt.gz",
    "A":   "{srr}_trimmomatic_adapter_fC.txt.gz",
    "P5":  "{srr}_trimmomatic_P5_fC.txt.gz",
    "P10": "{srr}_trimmomatic_P10_fC.txt.gz",
    "P20": "{srr}_trimmomatic_P20_fC.txt.gz",
    "P35": "{srr}_trimmomatic_P35_fC.txt.gz",
}
METHODS = list(METHOD_MAP.keys())


def load_count_file(path):
    """Load a featureCounts .gz file, returning a Series of gene → count."""
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    # Skip comment lines (start with #) and header
    data_lines = [l for l in lines if not l.startswith("#")]
    if not data_lines:
        return None
    header = data_lines[0].strip().split("\t")
    rows = [l.strip().split("\t") for l in data_lines[1:]]
    # featureCounts format: Geneid, Chr, Start, End, Strand, Length, Count
    genes = [r[0] for r in rows]
    counts = [int(r[-1]) for r in rows]
    return pd.Series(counts, index=genes, name=path.stem)


def build_count_matrix(project, method, srr_list):
    """Build count matrix (genes × samples) for a project and method."""
    proj_dir = COUNTS_DIR / project
    series_list = []
    valid_srrs = []

    for srr in srr_list:
        fname = METHOD_MAP[method].format(srr=srr)
        fpath = proj_dir / fname
        if not fpath.exists():
            continue
        s = load_count_file(fpath)
        if s is not None:
            s.name = srr
            series_list.append(s)
            valid_srrs.append(srr)

    if not series_list:
        return None, []
    mat = pd.concat(series_list, axis=1).fillna(0).astype(int)
    return mat, valid_srrs


def run_deseq2(count_matrix, conditions):
    """Run PyDESeq2 on a count matrix. Returns results DataFrame."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    metadata = pd.DataFrame({"condition": conditions}, index=count_matrix.columns)
    dds = DeseqDataSet(counts=count_matrix.T, metadata=metadata,
                       design="~condition")
    dds.deseq2()
    
    unique_conds = pd.Series(conditions).value_counts().index.tolist()
    if len(unique_conds) > 2:
        # For multi-class, contrast the two most frequent classes
        contrast = ["condition", unique_conds[0], unique_conds[1]]
        stat = DeseqStats(dds, contrast=contrast)
    else:
        stat = DeseqStats(dds)
        
    stat.summary()
    return stat.results_df


def run_gsea(de_results, gmt_path):
    """Run GSEA prerank using Wald statistics from DESeq2."""
    import gseapy as gp

    # Use Wald statistic for ranking
    rnk = de_results["stat"].dropna().sort_values(ascending=False)
    if len(rnk) < 100:
        return None

    res = gp.prerank(
        rnk=rnk,
        gene_sets=str(gmt_path),
        outdir=None,
        min_size=15,
        max_size=500,
        permutation_num=1000,
        seed=42,
        verbose=False,
    )
    return res.res2d


def discover_missing_projects():
    """Find projects that have sample sheets and counts but no DE/GSEA output."""
    missing = []
    for sheet in sorted(META_DIR.glob("PRJNA*.csv")):
        proj = sheet.stem
        proj_counts = COUNTS_DIR / proj
        proj_out = CONC_DIR / proj

        if not proj_counts.exists():
            continue
        # Check if already done
        if (proj_out / "de_U.tsv").exists() and (proj_out / "gsea_U.tsv").exists():
            continue
        missing.append(proj)
    return missing


def main():
    # Determine which projects to run
    if len(sys.argv) > 1:
        projects = sys.argv[1:]
    else:
        projects = discover_missing_projects()

    if not projects:
        print("All projects already have DE/GSEA results. Nothing to do.")
        return

    print(f"Will process {len(projects)} projects: {', '.join(projects)}")
    assert GMT_FILE.exists(), f"GMT file not found: {GMT_FILE}"

    for i, proj in enumerate(projects):
        print(f"\n{'='*60}")
        print(f"Project {i+1}/{len(projects)}: {proj}")
        print(f"{'='*60}")

        # Load sample sheet
        sheet_path = META_DIR / f"{proj}.csv"
        if not sheet_path.exists():
            print(f"  SKIP: no sample sheet at {sheet_path}")
            continue

        sheet = pd.read_csv(sheet_path)
        srr_list = sheet["Run"].tolist()
        conditions = sheet["condition"].tolist()

        if len(set(conditions)) < 2:
            print(f"  SKIP: only {len(set(conditions))} condition class(es)")
            continue

        print(f"  {len(srr_list)} samples, {len(set(conditions))} classes: "
              f"{pd.Series(conditions).value_counts().to_dict()}")

        out_dir = CONC_DIR / proj
        out_dir.mkdir(parents=True, exist_ok=True)

        for m in METHODS:
            print(f"    Method {m}...", end="", flush=True)
            try:
                count_matrix, valid_srrs = build_count_matrix(proj, m, srr_list)
                if count_matrix is None or len(valid_srrs) < 3:
                    print(f" not enough count files ({len(valid_srrs) if valid_srrs else 0})")
                    continue

                # Match conditions to valid SRRs
                srr_to_cond = dict(zip(srr_list, conditions))
                valid_conds = [srr_to_cond[s] for s in valid_srrs]

                if len(set(valid_conds)) < 2:
                    print(" only 1 condition in valid samples")
                    continue

                # DESeq2
                de_results = run_deseq2(count_matrix, valid_conds)
                de_results.to_csv(out_dir / f"de_{m}.tsv", sep="\t")

                # GSEA
                gsea_results = run_gsea(de_results, GMT_FILE)
                if gsea_results is not None:
                    gsea_results.to_csv(out_dir / f"gsea_{m}.tsv", sep="\t")
                    print(f" done ({len(de_results)} genes, "
                          f"{len(gsea_results)} pathways)")
                else:
                    print(f" done ({len(de_results)} genes, GSEA skipped)")

            except Exception as e:
                print(f" FAILED: {e}")

    print(f"\nFinished. Results in {CONC_DIR}/")


if __name__ == "__main__":
    main()
