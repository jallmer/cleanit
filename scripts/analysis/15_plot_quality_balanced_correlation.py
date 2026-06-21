#!/usr/bin/env python3
"""Plot balanced quality ranges against count-correlation loss by cleaning mode."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ANALYSIS = Path("/scratch/hpc-prf-omiks/ja/analysis")
PLOTS = ANALYSIS / "plots"
RANDOM_SEED = 17

PEARSON_MODES = {
    "Adapter": "untrmd_adptrTrmd_pear",
    "P5": "untrmd_P5Trmd_pear",
    "P10": "untrmd_P10Trmd_pear",
    "P20": "untrmd_P20Trmd_pear",
}


def load_inputs() -> pd.DataFrame:
    eval_df = pd.read_csv(ANALYSIS / "per_srr_eval.tsv", sep="\t")
    quality = pd.read_csv(ANALYSIS / "per_srr_quality.tsv", sep="\t")
    platform = pd.read_csv(ANALYSIS / "current_srr_platform_check.tsv", sep="\t")

    keep_quality = [
        "SRR_ID",
        "project_id",
        "Q_mean",
        "frac_below_q30",
        "tail_quality_decay",
        "read_length_mean",
        "duplication_rate",
    ]
    keep_platform = ["SRR_ID", "project_id", "library_strategy", "instrument_model"]
    df = eval_df.merge(quality[keep_quality], on=["SRR_ID", "project_id"], how="inner")
    df = df.merge(platform[keep_platform], on=["SRR_ID", "project_id"], how="left")

    numeric_cols = list(PEARSON_MODES.values()) + keep_quality[2:]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["Q_mean", *PEARSON_MODES.values()])
    return df


def assign_balanced_quality_bins(df: pd.DataFrame, bins: int = 4) -> tuple[pd.DataFrame, list[str]]:
    work = df.copy()
    work["_quality_bin_id"] = pd.qcut(work["Q_mean"], q=bins, labels=False, duplicates="drop")
    labels = []
    for bin_id in sorted(work["_quality_bin_id"].dropna().unique()):
        in_bin = work.loc[work["_quality_bin_id"] == bin_id, "Q_mean"]
        labels.append(f"Q{int(bin_id) + 1}: {in_bin.min():.1f}-{in_bin.max():.1f}")
    label_map = {bin_id: labels[i] for i, bin_id in enumerate(sorted(work["_quality_bin_id"].dropna().unique()))}
    work["quality_range"] = work["_quality_bin_id"].map(label_map)

    counts = work["quality_range"].value_counts()
    if counts.empty:
        return work.iloc[0:0].copy(), labels
    target_n = int(counts.min())
    sampled = (
        work.groupby("quality_range", group_keys=False, observed=True)
        .apply(lambda g: g.sample(n=target_n, random_state=RANDOM_SEED))
        .reset_index(drop=True)
    )
    return sampled, labels


def to_long(sampled: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mode, col in PEARSON_MODES.items():
        part = sampled[
            [
                "SRR_ID",
                "project_id",
                "library_strategy",
                "quality_range",
                "Q_mean",
                "read_length_mean",
                "duplication_rate",
                col,
            ]
        ].copy()
        part = part.rename(columns={col: "pearson"})
        part["cleaning_mode"] = mode
        part["correlation_loss"] = 1.0 - part["pearson"]
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def plot_subset(df: pd.DataFrame, subset_name: str, title_suffix: str) -> dict[str, int]:
    sampled, labels = assign_balanced_quality_bins(df)
    long = to_long(sampled)
    long.to_csv(ANALYSIS / f"quality_balanced_correlation_loss_{subset_name}.tsv", sep="\t", index=False)

    modes = list(PEARSON_MODES.keys())
    colors = ["#4C78A8", "#59A14F", "#F2CF5B", "#E15759"]
    x = np.arange(len(modes))
    width = 0.18

    fig, ax = plt.subplots(figsize=(11, 6))
    for i, label in enumerate(labels):
        vals = [
            long.loc[(long["quality_range"] == label) & (long["cleaning_mode"] == mode), "correlation_loss"].dropna()
            for mode in modes
        ]
        pos = x + (i - (len(labels) - 1) / 2) * width
        bp = ax.boxplot(
            vals,
            positions=pos,
            widths=width * 0.85,
            patch_artist=True,
            showfliers=False,
            manage_ticks=False,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(colors[i % len(colors)])
            patch.set_alpha(0.65)
        for median in bp["medians"]:
            median.set_color("#222222")
            median.set_linewidth(1.2)
        means = [v.mean() if len(v) else np.nan for v in vals]
        ax.plot(pos, means, color=colors[i % len(colors)], marker="o", linewidth=1.5, label=label)

    ax.set_xticks(x)
    ax.set_xticklabels(modes)
    ax.set_ylabel("Correlation loss vs untrimmed counts (1 - Pearson r)")
    ax.set_xlabel("Cleaning mode")
    ax.set_title(f"Balanced quality ranges vs count correlation loss ({title_suffix})")
    ax.set_yscale("symlog", linthresh=1e-5)
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(title="Q_mean range", fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOTS / f"quality_balanced_correlation_loss_{subset_name}.png", dpi=170)
    plt.close(fig)

    max_loss = sampled.copy()
    max_loss["max_correlation_loss_no_p35"] = 1.0 - max_loss[list(PEARSON_MODES.values())].max(axis=1)
    fig, ax = plt.subplots(figsize=(8, 5))
    grouped = [max_loss.loc[max_loss["quality_range"] == label, "max_correlation_loss_no_p35"].dropna() for label in labels]
    bp = ax.boxplot(grouped, labels=labels, patch_artist=True, showfliers=True)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
    ax.set_yscale("symlog", linthresh=1e-5)
    ax.set_ylabel("Worst correlation loss across Adapter/P5/P10/P20")
    ax.set_xlabel("Balanced Q_mean range")
    ax.set_title(f"Worst non-P35 count-correlation loss by quality range ({title_suffix})")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOTS / f"quality_balanced_max_correlation_loss_{subset_name}.png", dpi=170)
    plt.close(fig)

    return {
        "input_srrs_with_quality_and_counts": int(len(df)),
        "balanced_srrs": int(len(sampled)),
        "bins": int(len(labels)),
        "per_bin_srrs": int(len(sampled) / len(labels)) if labels else 0,
    }


def main() -> int:
    PLOTS.mkdir(parents=True, exist_ok=True)
    df = load_inputs()
    summary_rows = []
    rnaseq = df.loc[df["library_strategy"] == "RNA-Seq"].copy()
    odd_q3_cluster = {
        "SRR35798696",
        "SRR35798697",
        "SRR35798698",
        "SRR35798695",
        "SRR35798688",
        "SRR35798687",
        "SRR6334449",
    }
    rnaseq_filtered = rnaseq.loc[
        (~rnaseq["SRR_ID"].isin(odd_q3_cluster))
        & (rnaseq["read_length_mean"] >= 50)
    ].copy()

    exclusions = rnaseq.loc[~rnaseq.index.isin(rnaseq_filtered.index)].copy()
    exclusions["exclusion_reason"] = ""
    exclusions.loc[exclusions["SRR_ID"].isin(odd_q3_cluster), "exclusion_reason"] = "odd_q3_correlation_loss_cluster"
    exclusions.loc[exclusions["read_length_mean"] < 50, "exclusion_reason"] = (
        exclusions.loc[exclusions["read_length_mean"] < 50, "exclusion_reason"]
        .replace("", "short_read_length_lt_50bp")
        .mask(lambda x: x != "short_read_length_lt_50bp", lambda x: x + ";short_read_length_lt_50bp")
    )
    exclusions.to_csv(ANALYSIS / "quality_balanced_correlation_loss_excluded_srrs.tsv", sep="\t", index=False)

    for subset_name, subset_df, title in [
        ("all_illumina", df, "all Illumina strategies with available QC"),
        ("rnaseq", rnaseq, "RNA-seq only with available QC"),
        ("rnaseq_filtered", rnaseq_filtered, "RNA-seq only, odd samples excluded"),
    ]:
        stats = plot_subset(subset_df, subset_name, title)
        stats["subset"] = subset_name
        summary_rows.append(stats)

    pd.DataFrame(summary_rows)[["subset", "input_srrs_with_quality_and_counts", "bins", "per_bin_srrs", "balanced_srrs"]].to_csv(
        ANALYSIS / "quality_balanced_correlation_loss_summary.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(summary_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
