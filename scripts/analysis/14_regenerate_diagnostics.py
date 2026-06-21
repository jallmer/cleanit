#!/usr/bin/env python3
"""Regenerate post-summary diagnostic tables and plots from current analysis TSVs."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path("/scratch/hpc-prf-omiks/ja")
ANALYSIS = BASE / "analysis"
PLOTS = ANALYSIS / "plots"
FASTQC_BASE = BASE / "flattened_fastqc_raw"

JSD_COLS = [
    "untrmd_adptrTrmd_jsd",
    "untrmd_P5Trmd_jsd",
    "untrmd_P10Trmd_jsd",
    "untrmd_P20Trmd_jsd",
    "untrmd_P35Trmd_jsd",
    "P5Trmd_P35Trmd_jsd",
]
JSD_NO_P35 = [
    "untrmd_adptrTrmd_jsd",
    "untrmd_P5Trmd_jsd",
    "untrmd_P10Trmd_jsd",
    "untrmd_P20Trmd_jsd",
]
PEAR_COLS = [
    "untrmd_adptrTrmd_pear",
    "untrmd_P5Trmd_pear",
    "untrmd_P10Trmd_pear",
    "untrmd_P20Trmd_pear",
    "untrmd_P35Trmd_pear",
    "P5Trmd_P35Trmd_pear",
]
PEAR_NO_P35 = [
    "untrmd_adptrTrmd_pear",
    "untrmd_P5Trmd_pear",
    "untrmd_P10Trmd_pear",
    "untrmd_P20Trmd_pear",
]
LABELS_ALL = ["Adapter", "P5", "P10", "P20", "P35", "P5 vs P35"]
LABELS_NO_P35 = ["Adapter", "P5", "P10", "P20"]


def numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def save_boxplot(df: pd.DataFrame, cols: list[str], labels: list[str], title: str, ylabel: str, path: Path) -> None:
    data = [df[col].dropna().to_numpy() for col in cols]
    fig, ax = plt.subplots(figsize=(12, 6))
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=True, widths=0.6)
    colors = ["#4C78A8", "#59A14F", "#F2CF5B", "#E15759", "#B07AA1", "#9C755F"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_project_heatmap(df: pd.DataFrame, cols: list[str], labels: list[str], title: str, cmap: str, path: Path) -> None:
    grouped = df.groupby("project_id")[cols].mean(numeric_only=True)
    grouped = grouped.loc[grouped.index.sort_values()]
    fig, ax = plt.subplots(figsize=(10, max(5, len(grouped) * 0.32)))
    matrix = grouped.to_numpy(dtype=float)
    im = ax.imshow(matrix, aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_yticks(range(len(grouped.index)))
    ax.set_yticklabels(grouped.index, fontsize=7)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.75)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_profiles(df: pd.DataFrame, cols: list[str], labels: list[str], title: str, ylabel: str, path: Path) -> None:
    order = df[cols].max(axis=1).sort_values().index
    fig, ax = plt.subplots(figsize=(13, 6))
    x = np.arange(len(order))
    for col, label in zip(cols, labels):
        ax.plot(x, df.loc[order, col].to_numpy(dtype=float), linewidth=0.8, alpha=0.8, label=label)
    ax.set_title(title)
    ax.set_xlabel("SRRs sorted by maximum divergence")
    ax.set_ylabel(ylabel)
    ax.legend(ncol=min(len(labels), 4), fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def throughput_no_p35() -> None:
    detail_path = ANALYSIS / "throughput_detail.tsv"
    if not detail_path.exists():
        return
    df = pd.read_csv(detail_path, sep="\t")
    df = numeric(df, ["duration_sec", "cores_used", "mb_per_sec_per_core", "total_reads_M"])
    df["core_seconds"] = df["duration_sec"] * df["cores_used"]
    stage = df["stage"].fillna("").astype(str)
    mode = df["mode"].fillna("").astype(str)
    keep = ~(stage.str.contains("P35", case=False) | mode.str.fullmatch("P35", case=False))
    df = df.loc[keep].copy()
    df["mode_label"] = df["mode"].fillna("").replace("", "all")

    by_cat = df.groupby("category", dropna=False)["core_seconds"].sum().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    by_cat.plot(kind="bar", ax=ax, color="#4C78A8")
    ax.set_title("Core-seconds by stage, excluding P35")
    ax.set_ylabel("core-seconds")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOTS / "stage_core_seconds_no_p35.png", dpi=160)
    plt.close(fig)

    plot_df = df.dropna(subset=["mb_per_sec_per_core"])
    groups = list(plot_df.groupby(["category", "mode_label"], sort=True))
    labels = [f"{cat}\n{mode}" for (cat, mode), _ in groups]
    data = [g["mb_per_sec_per_core"].to_numpy(dtype=float) for _, g in groups]
    fig, ax = plt.subplots(figsize=(max(10, len(labels) * 0.55), 6))
    ax.boxplot(data, labels=labels, showfliers=False)
    ax.set_title("Throughput per core, excluding P35")
    ax.set_ylabel("MB/s/core")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOTS / "throughput_per_core_no_p35.png", dpi=160)
    plt.close(fig)

    scatter = df.dropna(subset=["total_reads_M", "core_seconds"])
    fig, ax = plt.subplots(figsize=(8, 6))
    for cat, g in scatter.groupby("category"):
        ax.scatter(g["total_reads_M"], g["core_seconds"], s=16, alpha=0.6, label=cat)
    ax.set_title("Core-seconds vs reads, excluding P35")
    ax.set_xlabel("Reads (millions)")
    ax.set_ylabel("core-seconds")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOTS / "core_seconds_vs_reads_no_p35.png", dpi=160)
    plt.close(fig)


def parse_bin_end(token: str) -> float:
    token = token.strip()
    if "-" in token:
        return float(token.split("-", 1)[1])
    return float(token)


def parse_fastqc_end_quality(path: Path) -> dict[str, float] | None:
    rows: list[list[str]] = []
    in_section = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">>Per base sequence quality"):
            in_section = True
            continue
        if in_section and line.startswith(">>END_MODULE"):
            break
        if in_section and line and not line.startswith("#"):
            rows.append(line.split("\t"))
    if not rows:
        return None
    last = rows[-1]
    tail = rows[-min(10, len(rows)) :]
    return {
        "last_bin_end": parse_bin_end(last[0]),
        "last_bin_q_mean": float(last[1]),
        "last_bin_q_median": float(last[2]),
        "end_q_mean": float(np.mean([float(r[1]) for r in tail])),
        "end_q_median": float(np.mean([float(r[2]) for r in tail])),
    }


def end_quality_table(outliers: pd.DataFrame) -> pd.DataFrame:
    wanted = set(outliers["SRR_ID"])
    rows = []
    for data_file in FASTQC_BASE.glob("*/*_fastqc/fastqc_data.txt"):
        fastqc_dir = data_file.parent.name
        srr = fastqc_dir.removesuffix("_fastqc")
        if srr.endswith("_1") or srr.endswith("_2"):
            srr = srr[:-2]
        if srr not in wanted:
            continue
        parsed = parse_fastqc_end_quality(data_file)
        if not parsed:
            continue
        row = {"project_id": data_file.parent.parent.name, "SRR_ID": srr}
        row.update(parsed)
        rows.append(row)
    endq = pd.DataFrame(rows)
    if endq.empty:
        return endq
    endq = endq.groupby(["project_id", "SRR_ID"], as_index=False).mean(numeric_only=True)
    return outliers.merge(endq, on=["project_id", "SRR_ID"], how="left")


def matched_controls(profile: pd.DataFrame, quality: pd.DataFrame) -> pd.DataFrame:
    qc_cols = [
        "sequence_depth",
        "Q_mean",
        "read_length_mean",
        "frac_below_q30",
        "tail_quality_decay",
        "adapter_rate",
        "duplication_rate",
        "gc_content",
    ]
    merged = profile.merge(quality[["SRR_ID", "project_id", *qc_cols]], on=["SRR_ID", "project_id"], how="inner")
    merged = numeric(merged, qc_cols)
    merged = merged.dropna(subset=qc_cols)
    if merged.empty:
        return merged
    features = pd.DataFrame(
        {
            "log_depth": np.log10(merged["sequence_depth"].clip(lower=1)),
            "Q_mean": merged["Q_mean"],
            "log_read_length": np.log10(merged["read_length_mean"].clip(lower=1)),
            "frac_below_q30": merged["frac_below_q30"],
            "tail_quality_decay": merged["tail_quality_decay"],
            "adapter_rate": merged["adapter_rate"],
            "duplication_rate": merged["duplication_rate"],
            "gc_content": merged["gc_content"],
        },
        index=merged.index,
    )
    z = (features - features.mean()) / features.std(ddof=0).replace(0, 1)
    bad_idx = merged.index[merged["max_no_p35_jsd"] >= 0.1].to_list()
    ctrl_idx = merged.index[merged["max_no_p35_jsd"] < 0.01].to_list()
    rows = []
    for idx in bad_idx:
        if not ctrl_idx:
            break
        dists = ((z.loc[ctrl_idx] - z.loc[idx]) ** 2).sum(axis=1) ** 0.5
        for rank, ctrl in enumerate(dists.sort_values().head(3).index, start=1):
            bad = merged.loc[idx]
            good = merged.loc[ctrl]
            rows.append(
                {
                    "bad_SRR_ID": bad["SRR_ID"],
                    "bad_project_id": bad["project_id"],
                    "control_rank": rank,
                    "control_SRR_ID": good["SRR_ID"],
                    "control_project_id": good["project_id"],
                    "qc_distance": dists.loc[ctrl],
                    "bad_max_no_p35_jsd": bad["max_no_p35_jsd"],
                    "control_max_no_p35_jsd": good["max_no_p35_jsd"],
                    "bad_read_length_mean": bad["read_length_mean"],
                    "control_read_length_mean": good["read_length_mean"],
                    "bad_duplication_rate": bad["duplication_rate"],
                    "control_duplication_rate": good["duplication_rate"],
                    "bad_Q_mean": bad["Q_mean"],
                    "control_Q_mean": good["Q_mean"],
                }
            )
    return pd.DataFrame(rows)


def save_qc_scatter(profile: pd.DataFrame, quality: pd.DataFrame, path: Path) -> None:
    merged = profile.merge(quality, on=["SRR_ID", "project_id"], how="inner")
    merged = numeric(merged, ["read_length_mean", "duplication_rate", "Q_mean", "max_no_p35_jsd"])
    merged = merged.dropna(subset=["read_length_mean", "duplication_rate", "max_no_p35_jsd"])
    fig, ax = plt.subplots(figsize=(8, 6))
    bad = merged["max_no_p35_jsd"] >= 0.1
    ax.scatter(merged.loc[~bad, "read_length_mean"], merged.loc[~bad, "duplication_rate"], s=18, alpha=0.45, color="#4C78A8", label="other")
    ax.scatter(merged.loc[bad, "read_length_mean"], merged.loc[bad, "duplication_rate"], s=30, alpha=0.9, color="#E15759", label="bad")
    ax.set_xlabel("Mean read length")
    ax.set_ylabel("Duplication rate")
    ax.set_title("Count-deviation outliers vs QC-matched population")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_qc_pca(profile: pd.DataFrame, quality: pd.DataFrame, path: Path) -> None:
    cols = ["sequence_depth", "Q_mean", "read_length_mean", "frac_below_q30", "tail_quality_decay", "adapter_rate", "duplication_rate", "gc_content"]
    merged = profile.merge(quality[["SRR_ID", "project_id", *cols]], on=["SRR_ID", "project_id"], how="inner")
    merged = numeric(merged, cols)
    merged = merged.dropna(subset=cols)
    if len(merged) < 3:
        return
    x = pd.DataFrame(
        {
            "log_depth": np.log10(merged["sequence_depth"].clip(lower=1)),
            "Q_mean": merged["Q_mean"],
            "log_read_length": np.log10(merged["read_length_mean"].clip(lower=1)),
            "frac_below_q30": merged["frac_below_q30"],
            "tail_quality_decay": merged["tail_quality_decay"],
            "adapter_rate": merged["adapter_rate"],
            "duplication_rate": merged["duplication_rate"],
            "gc_content": merged["gc_content"],
        }
    )
    z = ((x - x.mean()) / x.std(ddof=0).replace(0, 1)).to_numpy()
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    pcs = z @ vt[:2].T
    bad = merged["max_no_p35_jsd"].to_numpy() >= 0.1
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(pcs[~bad, 0], pcs[~bad, 1], s=18, alpha=0.45, color="#4C78A8", label="other")
    ax.scatter(pcs[bad, 0], pcs[bad, 1], s=30, alpha=0.9, color="#E15759", label="bad")
    ax.set_xlabel("QC PC1")
    ax.set_ylabel("QC PC2")
    ax.set_title("QC PCA of count-deviation outliers")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def count_eval_diagnostics() -> None:
    eval_path = ANALYSIS / "per_srr_eval.tsv"
    quality_path = ANALYSIS / "per_srr_quality.tsv"
    if not eval_path.exists():
        return
    df = pd.read_csv(eval_path, sep="\t")
    df = numeric(df, JSD_COLS + PEAR_COLS)

    save_boxplot(df, JSD_COLS, LABELS_ALL, "Count matrix JSD across trimming modes", "JSD", PLOTS / "count_eval_boxplot_jsd_all_srrs.png")
    save_boxplot(df, JSD_NO_P35, LABELS_NO_P35, "Count matrix JSD, excluding P35", "JSD", PLOTS / "count_eval_boxplot_jsd_no_p35.png")
    save_boxplot(df, PEAR_COLS, LABELS_ALL, "Count matrix Pearson across trimming modes", "Pearson r", PLOTS / "count_eval_boxplot_pearson_all_srrs.png")
    save_boxplot(df, PEAR_NO_P35, LABELS_NO_P35, "Count matrix Pearson, excluding P35", "Pearson r", PLOTS / "count_eval_boxplot_pearson_no_p35.png")
    save_project_heatmap(df, JSD_COLS, LABELS_ALL, "Mean JSD by project and mode", "magma", PLOTS / "count_eval_heatmap_jsd_by_project.png")
    save_project_heatmap(df, JSD_NO_P35, LABELS_NO_P35, "Mean JSD by project and mode, excluding P35", "magma", PLOTS / "count_eval_heatmap_jsd_by_project_no_p35.png")
    save_project_heatmap(df, PEAR_COLS, LABELS_ALL, "Mean Pearson by project and mode", "RdYlGn", PLOTS / "count_eval_heatmap_pearson_by_project.png")
    save_project_heatmap(df, PEAR_NO_P35, LABELS_NO_P35, "Mean Pearson by project and mode, excluding P35", "RdYlGn", PLOTS / "count_eval_heatmap_pearson_by_project_no_p35.png")
    save_profiles(df, JSD_COLS, LABELS_ALL, "SRR JSD profiles across trimming modes", "JSD", PLOTS / "count_eval_profile_jsd_by_srr.png")
    save_profiles(df, JSD_NO_P35, LABELS_NO_P35, "SRR JSD profiles, excluding P35", "JSD", PLOTS / "count_eval_profile_jsd_by_srr_no_p35.png")
    save_profiles(df, PEAR_COLS, LABELS_ALL, "SRR Pearson profiles across trimming modes", "Pearson r", PLOTS / "count_eval_profile_pearson_by_srr.png")
    save_profiles(df, PEAR_NO_P35, LABELS_NO_P35, "SRR Pearson profiles, excluding P35", "Pearson r", PLOTS / "count_eval_profile_pearson_by_srr_no_p35.png")

    profile = df[["SRR_ID", "project_id", *JSD_NO_P35]].copy()
    profile["max_no_p35_jsd"] = profile[JSD_NO_P35].max(axis=1)
    profile["worst_no_p35_comparison"] = profile[JSD_NO_P35].idxmax(axis=1)
    profile = profile.sort_values("max_no_p35_jsd", ascending=False)
    profile.to_csv(ANALYSIS / "count_eval_outlier_profile.tsv", sep="\t", index=False)

    top = profile.head(50).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, max(6, len(top) * 0.16)))
    ax.barh(top["SRR_ID"], top["max_no_p35_jsd"], color="#E15759")
    ax.set_xlabel("Max JSD excluding P35")
    ax.set_title("Top count-deviation outliers, excluding P35")
    fig.tight_layout()
    fig.savefig(PLOTS / "count_eval_top_outliers_jsd_no_p35.png", dpi=160)
    plt.close(fig)

    if quality_path.exists():
        quality = pd.read_csv(quality_path, sep="\t")
        quality = numeric(quality, ["sequence_depth", "Q_mean", "read_length_mean", "frac_below_q30", "tail_quality_decay", "adapter_rate", "duplication_rate", "gc_content"])
        qmerged = profile.merge(quality, on=["SRR_ID", "project_id"], how="inner")
        for xcol, outname in [
            ("adapter_rate", "count_eval_qc_scatter_untrimmed_vs_adapter_jsd.png"),
            ("read_length_mean", "count_eval_qc_scatter_untrimmed_vs_p20_jsd.png"),
        ]:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(qmerged[xcol], qmerged["max_no_p35_jsd"], s=18, alpha=0.55, color="#4C78A8")
            ax.set_xlabel(xcol)
            ax.set_ylabel("Max JSD excluding P35")
            ax.set_title(f"Max count deviation vs {xcol}")
            ax.grid(alpha=0.25)
            fig.tight_layout()
            fig.savefig(PLOTS / outname, dpi=160)
            plt.close(fig)

        endq = end_quality_table(profile)
        if not endq.empty:
            endq.to_csv(ANALYSIS / "count_eval_outlier_end_quality.tsv", sep="\t", index=False)

        controls = matched_controls(profile, quality)
        controls.to_csv(ANALYSIS / "count_eval_bad_vs_matched_controls.tsv", sep="\t", index=False)
        save_qc_scatter(profile, quality, PLOTS / "count_eval_bad_vs_matched_controls_qc_scatter.png")
        save_qc_pca(profile, quality, PLOTS / "count_eval_bad_vs_matched_controls_qc_pca.png")

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(profile["max_no_p35_jsd"].dropna(), bins=50, color="#4C78A8", alpha=0.8)
        ax.axvline(0.1, color="#E15759", linestyle="--", label="bad threshold")
        ax.axvline(0.01, color="#59A14F", linestyle="--", label="control threshold")
        ax.set_xlabel("Max JSD excluding P35")
        ax.set_ylabel("SRRs")
        ax.set_title("Distribution of non-P35 count deviations")
        ax.legend()
        fig.tight_layout()
        fig.savefig(PLOTS / "count_eval_bad_vs_matched_controls_jsd_hist.png", dpi=160)
        plt.close(fig)


def current_platform_summary() -> None:
    platform_path = ANALYSIS / "srr_platform_check.tsv"
    eval_path = ANALYSIS / "per_srr_eval.tsv"
    if not (platform_path.exists() and eval_path.exists()):
        return
    platform = pd.read_csv(platform_path, sep="\t")
    platform = platform.rename(
        columns={
            "run_accession": "SRR_ID",
            "study_accession": "project_id",
            "instrument_platform": "platform",
        }
    )
    current = pd.read_csv(eval_path, sep="\t", usecols=["SRR_ID", "project_id"])
    joined = current.merge(platform, on=["SRR_ID", "project_id"], how="left")
    summary = (
        joined.groupby(["platform", "instrument_model", "library_strategy"], dropna=False)
        .size()
        .reset_index(name="current_evaluated_srrs")
        .sort_values(["platform", "instrument_model", "library_strategy"])
    )
    summary.to_csv(ANALYSIS / "current_platform_strategy_summary.tsv", sep="\t", index=False)
    joined.to_csv(ANALYSIS / "current_srr_platform_check.tsv", sep="\t", index=False)


def main() -> int:
    PLOTS.mkdir(parents=True, exist_ok=True)
    throughput_no_p35()
    count_eval_diagnostics()
    current_platform_summary()
    print("Regenerated diagnostics in", ANALYSIS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
