#!/usr/bin/env python3
"""Plot tail-quality bins vs count-correlation loss within read-length buckets."""

from __future__ import annotations

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
RANDOM_SEED = 23
KNOWN_SUSPICIOUS_SRRS = {
    # Project-specific cluster with elevated loss across multiple cleaning modes.
    "SRR35798687",
    "SRR35798688",
    "SRR35798695",
    "SRR35798696",
    "SRR35798697",
    "SRR35798698",
    # Shorter-read Q3 sample with P20-specific behavior.
    "SRR6334449",
    # Very short RNA-seq-like libraries are not comparable to standard mRNA-seq.
    "SRR7058572",
    "SRR7058573",
    "SRR7058574",
    "SRR7058575",
    "SRR7058576",
    "SRR7058577",
    "SRR7058578",
    "SRR7058579",
    "SRR7058580",
    "SRR19842866",
    "SRR17165227",
}

PEARSON_MODES = {
    "Adapter": "untrmd_adptrTrmd_pear",
    "P5": "untrmd_P5Trmd_pear",
    "P10": "untrmd_P10Trmd_pear",
    "P20": "untrmd_P20Trmd_pear",
}


def parse_range_end(token: str) -> float:
    token = token.strip()
    if "-" in token:
        return float(token.split("-", 1)[1])
    return float(token)


def parse_fastqc_tail(path: Path) -> dict[str, float] | None:
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

    tail = rows[-min(10, len(rows)) :]
    first = rows[: min(10, len(rows))]
    return {
        "tail_q_mean": float(np.mean([float(r[1]) for r in tail])),
        "tail_q_median": float(np.mean([float(r[2]) for r in tail])),
        "last_bin_q_mean": float(rows[-1][1]),
        "last_bin_q_median": float(rows[-1][2]),
        "last_bin_end": parse_range_end(rows[-1][0]),
        "head_to_tail_q_drop": float(np.mean([float(r[1]) for r in first]) - np.mean([float(r[1]) for r in tail])),
    }


def build_tail_quality_table() -> pd.DataFrame:
    rows = []
    for data_file in FASTQC_BASE.glob("*/*_fastqc/fastqc_data.txt"):
        fastqc_name = data_file.parent.name.removesuffix("_fastqc")
        srr = fastqc_name[:-2] if fastqc_name.endswith(("_1", "_2")) else fastqc_name
        parsed = parse_fastqc_tail(data_file)
        if parsed is None:
            continue
        row = {"project_id": data_file.parent.parent.name, "SRR_ID": srr}
        row.update(parsed)
        rows.append(row)

    tail = pd.DataFrame(rows)
    if tail.empty:
        return tail
    return tail.groupby(["project_id", "SRR_ID"], as_index=False).mean(numeric_only=True)


def read_length_bucket(length: float) -> str:
    if length < 75:
        return "<75 bp"
    if length < 125:
        return "75-124 bp"
    if length < 175:
        return "125-174 bp"
    return ">=175 bp"


def load_analysis() -> pd.DataFrame:
    eval_df = pd.read_csv(ANALYSIS / "per_srr_eval.tsv", sep="\t")
    quality = pd.read_csv(ANALYSIS / "per_srr_quality.tsv", sep="\t")
    platform = pd.read_csv(ANALYSIS / "current_srr_platform_check.tsv", sep="\t")
    tail = build_tail_quality_table()
    tail.to_csv(ANALYSIS / "per_srr_tail_quality.tsv", sep="\t", index=False)

    cols = ["SRR_ID", "project_id", "read_length_mean", "duplication_rate", "Q_mean"]
    df = eval_df.merge(quality[cols], on=["SRR_ID", "project_id"], how="inner")
    df = df.merge(tail, on=["SRR_ID", "project_id"], how="inner")
    df = df.merge(platform[["SRR_ID", "project_id", "library_strategy"]], on=["SRR_ID", "project_id"], how="left")

    for col in [*PEARSON_MODES.values(), "read_length_mean", "duplication_rate", "Q_mean", "tail_q_mean", "tail_q_median", "last_bin_q_mean", "last_bin_q_median", "head_to_tail_q_drop"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=[*PEARSON_MODES.values(), "read_length_mean", "tail_q_median"])
    df = df.loc[(df["library_strategy"] == "RNA-Seq") & (df["read_length_mean"] >= 50)].copy()
    df["read_length_bucket"] = df["read_length_mean"].map(read_length_bucket)
    return df


def balanced_tail_bins(df: pd.DataFrame, q_col: str = "tail_q_median", bins: int = 4) -> tuple[pd.DataFrame, list[str]]:
    work = df.copy()
    work["_tail_bin_id"] = pd.qcut(work[q_col], q=bins, labels=False, duplicates="drop")
    ids = sorted(work["_tail_bin_id"].dropna().unique())
    labels = []
    for bin_id in ids:
        values = work.loc[work["_tail_bin_id"] == bin_id, q_col]
        labels.append(f"T{int(bin_id) + 1}: {values.min():.1f}-{values.max():.1f}")
    label_map = {bin_id: labels[i] for i, bin_id in enumerate(ids)}
    work["tail_quality_range"] = work["_tail_bin_id"].map(label_map)
    counts = work["tail_quality_range"].value_counts()
    if counts.empty:
        return work.iloc[0:0].copy(), labels
    target_n = int(counts.min())
    sampled = (
        work.groupby("tail_quality_range", group_keys=False, observed=True)
        .apply(lambda g: g.sample(n=target_n, random_state=RANDOM_SEED))
        .reset_index(drop=True)
    )
    return sampled, labels


def long_loss(sampled: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for mode, col in PEARSON_MODES.items():
        part = sampled[
            [
                "SRR_ID",
                "project_id",
                "read_length_bucket",
                "tail_quality_range",
                "tail_q_median",
                "last_bin_q_median",
                "head_to_tail_q_drop",
                "read_length_mean",
                "duplication_rate",
                col,
            ]
        ].copy()
        part = part.rename(columns={col: "pearson"})
        part["cleaning_mode"] = mode
        part["correlation_loss"] = 1.0 - part["pearson"]
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def plot_bucket(bucket: str, df: pd.DataFrame, suffix: str = "") -> dict[str, object] | None:
    sampled, labels = balanced_tail_bins(df)
    if sampled.empty or len(labels) < 2:
        return None

    safe_bucket = bucket.replace("<", "lt").replace(">=", "ge").replace(" ", "").replace("-", "_")
    long = long_loss(sampled)
    suffix_part = f"_{suffix}" if suffix else ""
    long.to_csv(ANALYSIS / f"tail_quality_correlation_loss_rnaseq_{safe_bucket}{suffix_part}.tsv", sep="\t", index=False)

    modes = list(PEARSON_MODES.keys())
    colors = ["#4C78A8", "#59A14F", "#F2CF5B", "#E15759"]
    x = np.arange(len(modes))
    width = min(0.18, 0.75 / len(labels))

    fig, ax = plt.subplots(figsize=(11, 6))
    for i, label in enumerate(labels):
        vals = [
            long.loc[(long["tail_quality_range"] == label) & (long["cleaning_mode"] == mode), "correlation_loss"].dropna()
            for mode in modes
        ]
        pos = x + (i - (len(labels) - 1) / 2) * width
        bp = ax.boxplot(vals, positions=pos, widths=width * 0.85, patch_artist=True, showfliers=False, manage_ticks=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(colors[i % len(colors)])
            patch.set_alpha(0.65)
        for median in bp["medians"]:
            median.set_color("#222222")
        means = [v.mean() if len(v) else np.nan for v in vals]
        ax.plot(pos, means, marker="o", linewidth=1.4, color=colors[i % len(colors)], label=label)

    ax.set_xticks(x)
    ax.set_xticklabels(modes)
    ax.set_yscale("symlog", linthresh=1e-5)
    ax.set_ylabel("Correlation loss vs untrimmed counts (1 - Pearson r)")
    ax.set_xlabel("Cleaning mode")
    title_suffix = " excluding suspicious SRRs" if suffix else ""
    ax.set_title(f"Tail quality vs count-correlation loss, RNA-seq {bucket}{title_suffix}")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(title="Tail median quality", fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOTS / f"tail_quality_correlation_loss_rnaseq_{safe_bucket}{suffix_part}.png", dpi=170)
    plt.close(fig)

    wide = sampled.copy()
    wide["max_correlation_loss"] = 1.0 - wide[list(PEARSON_MODES.values())].max(axis=1)
    fig, ax = plt.subplots(figsize=(8, 5))
    grouped = [wide.loc[wide["tail_quality_range"] == label, "max_correlation_loss"].dropna() for label in labels]
    bp = ax.boxplot(grouped, labels=labels, patch_artist=True, showfliers=True)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
    ax.set_yscale("symlog", linthresh=1e-5)
    ax.set_ylabel("Worst correlation loss across Adapter/P5/P10/P20")
    ax.set_xlabel("Balanced tail-quality range")
    ax.set_title(f"Worst non-P35 loss by tail quality, RNA-seq {bucket}{title_suffix}")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOTS / f"tail_quality_max_correlation_loss_rnaseq_{safe_bucket}{suffix_part}.png", dpi=170)
    plt.close(fig)

    return {
        "subset": suffix or "all_eligible",
        "read_length_bucket": bucket,
        "eligible_srrs": int(len(df)),
        "balanced_srrs": int(len(sampled)),
        "tail_quality_bins": int(len(labels)),
        "per_bin_srrs": int(len(sampled) / len(labels)),
        "tail_quality_ranges": "; ".join(labels),
    }


def main() -> int:
    PLOTS.mkdir(parents=True, exist_ok=True)
    df = load_analysis()
    rows = []
    nonsuspicious = df.loc[~df["SRR_ID"].isin(KNOWN_SUSPICIOUS_SRRS)].copy()
    excluded = df.loc[df["SRR_ID"].isin(KNOWN_SUSPICIOUS_SRRS)].copy()
    excluded["exclusion_reason"] = "known_suspicious_srr"
    excluded.to_csv(ANALYSIS / "tail_quality_correlation_loss_excluded_srrs.tsv", sep="\t", index=False)

    for subset_name, subset_df in [("all_eligible", df), ("nonsuspicious", nonsuspicious)]:
        suffix = "" if subset_name == "all_eligible" else subset_name
        for bucket in ["<75 bp", "75-124 bp", "125-174 bp", ">=175 bp"]:
            bucket_df = subset_df.loc[subset_df["read_length_bucket"] == bucket].copy()
            if len(bucket_df) < 12:
                continue
            result = plot_bucket(bucket, bucket_df, suffix=suffix)
            if result:
                rows.append(result)

    summary = pd.DataFrame(rows)
    summary.to_csv(ANALYSIS / "tail_quality_correlation_loss_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
