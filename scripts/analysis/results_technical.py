#!/usr/bin/env python
# coding: utf-8

# # Technical Results Notebook
# 
# Methodological details are kept separately in `methodology.ipynb`
# 
# This notebook contains the result narrative and the code used to aggregate the result tables and plots. 
# 
# Tables and figures are generated from the analysis-ready technical outputs, with methodological cohort construction documented in `methodology.ipynb`. Some operational details were stored in text summaries rather than retaining SAM/BAM or FASTQ files, which would have been too costly.
# 
# While processing more than 1000 SRR files, many questions arose that were not included originally. Therefore, some of the analyses here are on partial data such as timing analyses. However, depending on the question, the amount of results available should still be representable enough for the side-questions that we asked along the way.

# In[1]:


from pathlib import Path
import sqlite3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ANALYSIS = Path.cwd()
if not (ANALYSIS / "technical").exists():
    ANALYSIS = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis")
BASE = ANALYSIS.parent
TECH = ANALYSIS / "technical"
DB = Path("/pc2/users/o/omiks001/srr_queue.db")

pd.set_option("display.max_columns", 100)
pd.set_option("display.width", 180)
plt.rcParams["figure.dpi"] = 120

SUSPICIOUS_SRRS = set("""
SRR35798687 SRR35798688 SRR35798695 SRR35798696 SRR35798697 SRR35798698
SRR6334449 SRR7058572 SRR7058573 SRR7058574 SRR7058575 SRR7058576
SRR7058577 SRR7058578 SRR7058579 SRR7058580 SRR19842866 SRR17165227
""".split())
try:
    from IPython.display import display
except ImportError:
    display = print


# ## Executive Interpretation
# 
# This technical analysis asks whether read cleaning changes downstream count profiles. Illumina RNA-seq samples were used in this assessment although not all represent standard RNA-seq analysis type processing. The answer is that cleaning does not materially change the results for adapter-only, P5, P10, and P20. These remain well aligned with the untrimmed counts. That means the technical effect of routine cleaning is small, so the biological DESeq2/GSEA follow-up is expected to be at least as stable unless a project has unusual library behavior (tbk).
# 
# Two exceptions are useful rather than contradictory. P35 is deliberately aggressive and shows the expected breakdown in read retention and count-profile similarity. It was a stress-test, not a routine cleaning choice. The non-P35 points with very high correlation loss are concentrated in `PRJNA1108066`, a miRNA-seq project that was we 'forced' through an mRNA/exon counting workflow; those points reflect a library-strategy mismatch and unstable correlations against nearly empty untrimmed count vectors. With shorter reads, random alignment became predominant.
# 
# Endpoint-download records should be interpreted cautiously. In the recorded operational attempts, ENA HTTP had the highest success rate and highest median throughput. This does not prove that Aspera is generally worse, because the endpoint attempts were sparse and unbalanced, but our own DB records do not support the older assumption that Aspera was faster and more stable in this workflow.
# 

# In[ ]:





# ## 1. Download Endpoint and SRA Conversion Results
# 
# The SRR download endpoint reliability and download speed are reconstructed from the database that steered the processing on the HPC. SRA conversion is a separate local processing step, so its throughput is taken from separate records not stored in the DB. This separates network endpoint performance from post-download decompression/conversion computational costs.

# In[2]:


con = sqlite3.connect(str(DB))
try:
    fetch = pd.read_sql_query("SELECT * FROM fetch_attempts", con, parse_dates=["started_at", "finished_at"])
finally:
    con.close()
fetch["duration_sec"] = (fetch["finished_at"] - fetch["started_at"]).dt.total_seconds()
fetch["mb_downloaded"] = pd.to_numeric(fetch["bytes_downloaded"], errors="coerce") / 1e6
fetch["mb_per_sec"] = fetch["mb_downloaded"] / fetch["duration_sec"]
endpoint_labels = {"aws":"AWS download", "ena_aspera":"ENA Aspera download", "ena_http":"ENA HTTP download", "gcp":"GCP download", "ncbi":"NCBI SRA download"}
data_types = {"aws":"SRA", "ena_aspera":"FASTQ.gz", "ena_http":"FASTQ.gz", "gcp":"SRA", "ncbi":"SRA"}
rows = []
for endpoint, grp in fetch.groupby("endpoint"):
    success = grp.loc[grp["status"].eq("success")].copy()
    attempts = len(grp)
    rows.append({
        "front_end_operation": endpoint_labels.get(endpoint, endpoint),
        "data_type": data_types.get(endpoint, "unknown"),
        "successes": len(success),
        "attempts": attempts,
        "failed": int(grp["status"].eq("failed").sum()),
        "cancelled": int(grp["status"].astype(str).str.contains("cancel", case=False, na=False).sum()),
        "success_pct": 100 * len(success) / attempts if attempts else np.nan,
        "median_mb_s": success["mb_per_sec"].replace([np.inf, -np.inf], np.nan).median(),
    })
endpoint_summary = pd.DataFrame(rows).sort_values("front_end_operation")
throughput = pd.read_csv(str(TECH / "throughput_detail.tsv"), sep="	")
sra_conversion = throughput.loc[(throughput["category"].eq("sra_conversion")) & (throughput["duration_sec"].ge(60))].copy()
endpoint_summary.loc[len(endpoint_summary)] = {
    "front_end_operation": "SRA conversion to .fastq.gz",
    "data_type": "SRA to FASTQ.gz",
    "successes": np.nan,
    "attempts": np.nan,
    "failed": np.nan,
    "cancelled": np.nan,
    "success_pct": np.nan,
    "median_mb_s": sra_conversion["mb_per_sec"].median(),
}
endpoint_summary


# In[3]:


# Endpoint diagnostics: the database contains only recorded fetch attempts,
# not a complete audit trail for all SRRs processed in the project.
fetch_success = fetch.loc[
    fetch["status"].eq("success")
    & fetch["duration_sec"].gt(0)
    & fetch["mb_downloaded"].gt(0)
].copy()

endpoint_distribution = (
    fetch_success.groupby("endpoint")
    .agg(
        speed_n=("srr_id", "size"),
        median_mb_s=("mb_per_sec", "median"),
        q1_mb_s=("mb_per_sec", lambda s: s.quantile(0.25)),
        q3_mb_s=("mb_per_sec", lambda s: s.quantile(0.75)),
        min_mb_s=("mb_per_sec", "min"),
        max_mb_s=("mb_per_sec", "max"),
        median_downloaded_gb=("mb_downloaded", lambda s: s.median() / 1000),
    )
    .reset_index()
    .sort_values("endpoint")
)

# Paired comparisons are sparse because endpoints were fallback attempts rather than
# a balanced benchmark where every SRR was downloaded from every endpoint.
best_success = fetch_success.sort_values("mb_per_sec", ascending=False).drop_duplicates(["srr_id", "endpoint"])
endpoint_wide = best_success.pivot(index="srr_id", columns="endpoint", values="mb_per_sec")
paired_rows = []
for a, b in [("ena_http", "ena_aspera"), ("ena_http", "aws"), ("ena_http", "gcp"), ("ena_http", "ncbi")]:
    if a not in endpoint_wide or b not in endpoint_wide:
        continue
    pair = endpoint_wide[[a, b]].dropna()
    if pair.empty:
        paired_rows.append({"comparison": f"{a} vs {b}", "paired_srrs": 0})
        continue
    diff = pair[a] - pair[b]
    paired_rows.append({
        "comparison": f"{a} vs {b}",
        "paired_srrs": len(pair),
        f"{a}_wins": int((diff > 0).sum()),
        f"{b}_wins": int((diff < 0).sum()),
        "median_mb_s_difference": diff.median(),
    })
paired_endpoint_comparison = pd.DataFrame(paired_rows)

endpoint_errors = (
    fetch.loc[fetch["status"].ne("success")]
    .groupby(["endpoint", "status", "error_msg"], dropna=False)
    .size()
    .reset_index(name="n")
    .sort_values(["endpoint", "status", "error_msg"])
)

display(endpoint_distribution)
display(paired_endpoint_comparison)
display(endpoint_errors)


# ### Endpoint Downloads
# 
# The recorded endpoint data do not support the earlier working assumption that Aspera was faster or more stable in this workflow. Among the recorded attempts, ENA HTTP had both the highest success rate and the highest median throughput. ENA Aspera had fewer successful transfers, more failed/cancelled attempts, and a lower median throughput than ENA HTTP.
# 
# These endpoint data are sparse and unbalanced because the pipeline used endpoint fallback rather than a controlled benchmark design. Therefore, in our operational records ENA HTTP was the most reliable and fastest endpoint overall, while the available paired data are too limited for a general network benchmark claim.
# 

# SRA files save some disk space relative to converted FASTQ.gz, but that saving comes with extra conversion cost. The following paired-size summary uses only SRRs where both SRA and converted FASTQ.gz sizes are available in the database.

# In[4]:


con = sqlite3.connect(str(DB))
try:
    queue = pd.read_sql_query("SELECT project_id, srr_id, sra_size_bytes, fastq_size_bytes, size_gb FROM srr_queue", con)
finally:
    con.close()
paired_sizes = queue.dropna(subset=["sra_size_bytes", "fastq_size_bytes"]).copy()
paired_sizes = paired_sizes.loc[(paired_sizes["sra_size_bytes"] > 0) & (paired_sizes["fastq_size_bytes"] > 0)].copy()
paired_sizes["sra_gb"] = paired_sizes["sra_size_bytes"] / 1e9
paired_sizes["fastq_gb"] = paired_sizes["fastq_size_bytes"] / 1e9
paired_sizes["sra_to_fastq_ratio"] = paired_sizes["sra_size_bytes"] / paired_sizes["fastq_size_bytes"]
paired_sizes[["sra_gb", "fastq_gb", "sra_to_fastq_ratio"]].agg(["count", "median", "mean"])
summary = pd.Series({
    "n_runs_with_both_sizes": len(paired_sizes),
    "total_sra_tb": paired_sizes["sra_size_bytes"].sum() / 1e12,
    "total_fastq_tb": paired_sizes["fastq_size_bytes"].sum() / 1e12,
    "fastq_minus_sra_tb": (
        paired_sizes["fastq_size_bytes"].sum()
        - paired_sizes["sra_size_bytes"].sum()
    ) / 1e12,
    "sra_as_percent_of_fastq": (
        paired_sizes["sra_size_bytes"].sum()
        / paired_sizes["fastq_size_bytes"].sum()
        * 100
    ),
    "median_run_sra_as_percent_of_fastq": (
        paired_sizes["sra_to_fastq_ratio"].median() * 100
    ),
})
summary


# ### Download and SRA Conversion
# 
# The endpoint table separates network retrieval from local SRA conversion. Direct FASTQ.gz retrieval avoids the SRA conversion step, whereas SRA-producing endpoints require local conversion/decompression before the same processing workflow can begin.
# 
# From a computational-cost perspective, SRA was less attractive than direct FASTQ.gz retrieval for this analysis. In the paired-size subset, SRA saved a minuite amount of disk space relative to converted FASTQ.gz, but the saving came with an additional conversion step, extra scheduling complexity, and substantial CPU use. At the observed conversion throughput, SRA conversion is not a trivial pre-processing detail; it is a meaningful compute cost that has to be accounted for when scaling to hundreds or thousands of SRRs. The compute cost equals that of the bowtie2 alignment step.
# 
# The practical result is that direct FASTQ.gz availability is preferable for this kind of large technical benchmark. SRA is useful as an archival format, but if SRA-side services also exposed ready-to-use `.fastq.gz` files consistently, workflows like this would spend less cluster time on conversion and more time on the actual analysis. Considering the very small space saving we observed, SRA should be avoided.
# 

# ## 2. Processing Cost
# 
# Processing cost is recorded on a per step basis, e.g., for FastQC analysis. The time and the amount of cores used is recorded for several hundred of the analyses done in this project. The six different trimming - alignment - counting branches are however not comparable among the branches becahse different settings were used and for p35 this becomes especially obvious. Among SRR files the comparison of the same branches is of course possible. 
# 
# 

# In[5]:


throughput = pd.read_csv(str(TECH / "throughput_detail.tsv"), sep="	")
for col in ["duration_sec", "cores_used", "total_reads_M", "fastq_mb"]:
    throughput[col] = pd.to_numeric(throughput[col], errors="coerce")
throughput["core_seconds"] = throughput["duration_sec"] * throughput["cores_used"]
throughput["cpu_hours"] = throughput["core_seconds"] / 3600

align_modes = ["untrimmed", "adapter_only", "P5", "P10", "P20", "P35"]
stage_order = ["fastqc", "sra_conversion"] + [f"{prefix}_{mode}" for mode in align_modes for prefix in ["trim", "align", "count"]]
# `trim_untrimmed` does not exist; keep only stages that occur in the timing table.
stage_order = [stage for stage in stage_order if throughput["stage"].eq(stage).any()]

stage_cost = (
    throughput.loc[throughput["stage"].isin(stage_order)]
    .groupby(["category", "stage", "mode"], dropna=False, as_index=False)
    .agg(
        records=("core_seconds", "size"),
        median_core_seconds=("core_seconds", "median"),
        q1_core_seconds=("core_seconds", lambda s: s.quantile(0.25)),
        q3_core_seconds=("core_seconds", lambda s: s.quantile(0.75)),
        median_cpu_hours=("cpu_hours", "median"),
    )
)
stage_cost["stage"] = pd.Categorical(stage_cost["stage"], categories=stage_order, ordered=True)
stage_cost = stage_cost.sort_values("stage")
stage_cost


# In[6]:


from matplotlib.patches import Patch

plot_data = []
labels = []
for stage in stage_order:
    vals = throughput.loc[throughput["stage"].eq(stage), "core_seconds"].dropna()
    vals = vals.loc[vals.gt(0)]
    if len(vals):
        plot_data.append(vals.to_numpy())
        labels.append(stage.replace("_", "\n"))

print(f"Data points used for plot: {sum(len(v) for v in plot_data)}")
fig, ax = plt.subplots(figsize=(20, 8))
bp = ax.boxplot(plot_data, labels=labels, patch_artist=True, notch=True,
                widths=0.6, flierprops=dict(markersize=3, alpha=0.4))
cat_colors = {"fastqc": "#42A5F5", "sra": "#78909C", "trim": "#66BB6A", "align": "#FFA726", "count": "#AB47BC"}
for i, label_text in enumerate(labels):
    flat = label_text.replace("\n", "_")
    if flat.startswith("trim"):
        color = cat_colors["trim"]
    elif flat.startswith("align"):
        color = cat_colors["align"]
    elif flat.startswith("count"):
        color = cat_colors["count"]
    elif flat.startswith("fastqc"):
        color = cat_colors["fastqc"]
    else:
        color = cat_colors["sra"]
    bp["boxes"][i].set_facecolor(color)
    bp["boxes"][i].set_alpha(0.7)

ax.set_ylabel("Core-seconds (duration x allocated cores)")
ax.set_title("Computational cost per pipeline stage")
ax.tick_params(axis="x", rotation=45, labelsize=8)
ax.grid(axis="y", alpha=0.3)
ax.set_yscale("log")
ax.legend(handles=[
    Patch(facecolor="#42A5F5", alpha=0.7, label="FastQC (8 cores)"),
    Patch(facecolor="#78909C", alpha=0.7, label="SRA conversion (36 cores)"),
    Patch(facecolor="#66BB6A", alpha=0.7, label="Trimming (36 cores)"),
    Patch(facecolor="#FFA726", alpha=0.7, label="Alignment (36 cores)"),
    Patch(facecolor="#AB47BC", alpha=0.7, label="Counting (36 cores)"),
], loc="upper right")
fig.tight_layout()
plt.show()


# ### Processing Cost
# 
# The stage-cost plot should be read as the cost of the benchmark matrix, not as the cost of a single normal RNA-seq analysis. A normal analysis would run FastQC, one trimming choice, one alignment, and one counting step. Here, the repeated cleaning/alignment/counting branches were run to test whether cleaning changes the count profile, so the computational cost is intentionally multiplied by the experimental design.
# 
# Translated back to a standard one-pass RNA-seq workflow, FastQC plus one trimming mode is not negligible. One Bowtie2 alignment costs roughly the same order of compute as FastQC plus one trimming pass, while counting is neglegible by comparison. This matters for interpretation: cleaning itself is not free and together with quality check adds a significant amount of compute. When running other aligners such as STAR, this could be much more cost for QC and trimming than for alignment.
# 
# Considering the SRA conversion, it added as much computational cost as the bowtie2 alignment. 
# 

# ## 3. Count-profile Stability Results
# 
# For each SRR and cleaning mode, the cleaned count vector is compared against the untrimmed count vector. Pearson correlation captures linear count-vector similarity, while Jensen-Shannon divergence captures distributional change.

# In[7]:


eval_df = pd.read_csv(str(TECH / "per_srr_eval.tsv"), sep="	")
mode_map = {
    "Adapter": ("untrmd_adptrTrmd_pear", "untrmd_adptrTrmd_jsd"),
    "P5": ("untrmd_P5Trmd_pear", "untrmd_P5Trmd_jsd"),
    "P10": ("untrmd_P10Trmd_pear", "untrmd_P10Trmd_jsd"),
    "P20": ("untrmd_P20Trmd_pear", "untrmd_P20Trmd_jsd"),
    "P35": ("untrmd_P35Trmd_pear", "untrmd_P35Trmd_jsd"),
}
count_stability = pd.DataFrame([
    {
        "mode": mode,
        "srrs_with_values": len(eval_df[[pear_col, jsd_col]].dropna()),
        "median_pearson_r": eval_df[pear_col].median(skipna=True),
        "median_jsd": eval_df[jsd_col].median(skipna=True),
    }
    for mode, (pear_col, jsd_col) in mode_map.items()
])
count_stability


# In[8]:


fig, ax = plt.subplots(figsize=(7.5, 4.5))
print(f"Data points used for plot: {len(eval_df.dropna(subset=[col for col, _ in mode_map.values()], how='all'))}")
ax.boxplot([eval_df[pear_col].dropna() for pear_col, _ in mode_map.values()], labels=list(mode_map.keys()), showfliers=False)
ax.set_ylabel("Pearson r vs untrimmed counts")
ax.set_title("Count-profile stability by cleaning mode")
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
plt.show()


# ### Interpretation: Count-profile Stability
# 
# The count-profile comparison is the core technical result. Across the Illumina assessment set, adapter-only, P5, P10, and P20 cleaning remain well aligned with the untrimmed count vectors for most SRRs. This supports the central view from the executive summary: routine cleaning does not materially change the count profiles in the datasets where the mRNA-style counting workflow is appropriate.
# 
# This statement is intentionally about routine cleaning and count-profile stability, not about every dataset being an ideal standard RNA-seq experiment. Some projects enter the broad transcriptomic filter but behave differently because the library type or read structure is not well matched to exon-level mRNA counting. Those cases are handled as interpretation caveats rather than as evidence that ordinary cleaning is unstable.
# 
# P35 is deliberately aggressive and behaves as the stress-test case where cleaning can damage the count profile. It should be retained as a stress test but excluded from the main practical cleaning recommendation.
# 

# ## 4. Read-retention Results
# 
# Read retention is evaluated from all available mode-resolved Trimmomatic summaries. The original logs in the read-only `fb/omiks_project/results` tree provide most of the data, and the local `flattened_trimmomatic_stats` folder contributes additional SRR/mode combinations that are not present in the `fb` tree. The sources are additive after de-duplication by `project_id + srr_id + mode`.
# 
# For duplicate SRR/mode keys, the original `fb` log-derived value is preferred because it is parsed directly from the Trimmomatic log. Flattened-only keys are then added. The combined table is cached as `technical/trimmomatic_retention_combined.tsv`.
# 

# In[9]:


trim = pd.read_csv(str(TECH / "trimmomatic_detail.tsv"), sep="\t")
# Use the consolidated trimmomatic stats natively!
if "SRR_ID" in trim.columns:
    trim = trim.rename(columns={"SRR_ID": "srr_id"})
trim["surviving_pct"] = pd.to_numeric(trim["surviving_pct"], errors="coerce")
trim_modes = ["adapter_only", "P5", "P10", "P20", "P35"]
wide_trim = trim.loc[trim["mode"].isin(trim_modes)].pivot_table(
    index=["project_id", "srr_id"], columns="mode", values="surviving_pct", aggfunc="mean"
)
paired_trim_practical = wide_trim.dropna(subset=["adapter_only", "P5", "P10", "P20"])
paired_trim_all_modes = wide_trim.dropna(subset=trim_modes)

retention_coverage = pd.DataFrame({
    "source_or_subset": [
        "flattened trimmomatic stats (all projects)",
        "combined mode-resolved table, any mode",
        "combined with adapter/P5/P10/P20",
        "combined with adapter/P5/P10/P20/P35",
    ],
    "count": [
        len(trim),
        wide_trim.shape[0],
        paired_trim_practical.shape[0],
        paired_trim_all_modes.shape[0],
    ],
})

read_retention = pd.DataFrame([
    {
        "mode": mode,
        "srrs_with_mode": int(wide_trim[mode].notna().sum()),
        "paired_all_modes_srrs": len(paired_trim_all_modes),
        "median_surviving_pct": wide_trim[mode].median(),
        "mean_surviving_pct": wide_trim[mode].mean(),
        "q1_surviving_pct": wide_trim[mode].quantile(0.25),
        "q3_surviving_pct": wide_trim[mode].quantile(0.75),
        "min_surviving_pct": wide_trim[mode].min(),
        "max_surviving_pct": wide_trim[mode].max(),
    }
    for mode in trim_modes
])

display(retention_coverage)
display(read_retention)


# In[10]:


fig, ax = plt.subplots(figsize=(7.5, 4.5))
ax.boxplot(
    [paired_trim_all_modes[mode].dropna() for mode in trim_modes],
    labels=["Adapter", "P5", "P10", "P20", "P35"],
    showfliers=False,
)
ax.set_ylabel("Surviving reads (%)")
ax.set_title(f"Paired read retention across trimming modes (n={len(paired_trim_all_modes)})")
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
plt.show()


# ### Interpretation: Read Retention
# 
# The read-retention sources are additive. The original `fb` Trimmomatic logs provide most mode-resolved retention records, and the local flattened stats add SRR/mode combinations that are not present in the `fb` tree. After de-duplication, the combined table contains 1167 SRRs with at least one trimming mode and 1027 SRRs with Adapter, P5, P10, P20, and P35 all present.
# 
# For duplicate SRR/mode keys, the `fb` log-derived value is preferred because it is parsed directly from the Trimmomatic log. This matters because a few overlapping flattened adapter-only rows are zero while the original logs show normal near-100% retention, suggesting that those flattened rows can include failed or incomplete records.
# 
# With the combined mode-resolved table, Adapter/P5/P10/P20 retain nearly all reads for most SRRs. P35 remains the clear aggressive stress-test condition and is not comparable to routine cleaning modes.
# 

# ## 5. Integrated Quality Sensitivity Results
# 
# The integrated quality view combines `srr_quality_scores.tsv` with count-profile metrics. Suspicious SRRs are excluded. The terminal mean quality metric prefers flattened FastQC raw output when available and falls back to MultiQC per-base quality otherwise.

# In[11]:


quality = pd.read_csv(str(TECH / "srr_quality_scores.tsv"), sep="	")
merged = eval_df.merge(quality, on="SRR_ID", how="inner", suffixes=("", "_quality"))
pear_cols = {
    "Adapter": "untrmd_adptrTrmd_pear",
    "P5": "untrmd_P5Trmd_pear",
    "P10": "untrmd_P10Trmd_pear",
    "P20": "untrmd_P20Trmd_pear",
    "P35": "untrmd_P35Trmd_pear",
}
for col in list(pear_cols.values()) + ["fastqc_raw_tail10_mean_q", "multiqc_per_base_tail10_mean_q", "fastqc_raw_read_length_end", "multiqc_per_base_read_length_end"]:
    merged[col] = pd.to_numeric(merged[col], errors="coerce")
merged["terminal_mean_q"] = merged["fastqc_raw_tail10_mean_q"].combine_first(merged["multiqc_per_base_tail10_mean_q"])
merged["read_length_end"] = merged["fastqc_raw_read_length_end"].combine_first(merged["multiqc_per_base_read_length_end"])
merged = merged.loc[~merged["SRR_ID"].isin(SUSPICIOUS_SRRS)].copy()
merged = merged.dropna(subset=list(pear_cols.values()) + ["terminal_mean_q", "read_length_end"])

def read_length_bucket(x):
    if x < 75:
        return "<75 bp"
    if x < 125:
        return "75-124 bp"
    if x < 175:
        return "125-174 bp"
    return ">=175 bp"

merged["read_length_bucket"] = merged["read_length_end"].map(read_length_bucket)
parts = []
for mode, col in pear_cols.items():
    part = merged[["SRR_ID", "project_id", "terminal_mean_q", "read_length_end", "read_length_bucket", col]].rename(columns={col: "pearson"}).copy()
    part["cleaning_mode"] = mode
    part["correlation_loss"] = 1.0 - part["pearson"]
    parts.append(part)
integrated_tail = pd.concat(parts, ignore_index=True)
quality_summary = merged.groupby("read_length_bucket")["terminal_mean_q"].agg(["count", "min", "median", "max"]).reset_index()
print(f"SRRs: {len(merged)}")
print(f"Mode-specific points: {len(integrated_tail)}")
quality_summary


# In[12]:


colors = {"<75 bp": "#4C78A8", "75-124 bp": "#59A14F", "125-174 bp": "#F28E2B", ">=175 bp": "#E15759"}
print(f"Data points used for plot: {len(integrated_tail)}")
fig, axes = plt.subplots(1, 5, figsize=(16, 4.8), sharex=True, sharey=True)
for ax, mode in zip(axes, pear_cols.keys()):
    sub = integrated_tail.loc[integrated_tail["cleaning_mode"].eq(mode)]
    for bucket, color in colors.items():
        ss = sub.loc[sub["read_length_bucket"].eq(bucket)]
        ax.scatter(ss["terminal_mean_q"], ss["correlation_loss"], s=13, alpha=0.45, color=color, label=bucket, edgecolors="none")
    ax.set_title(mode)
    ax.set_yscale("symlog", linthresh=1e-5)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.22)
    ax.set_xlabel("Terminal mean Q")
axes[0].set_ylabel("Correlation loss vs untrimmed (1 - Pearson r)")
handles, labels = axes[-1].get_legend_handles_labels()
fig.legend(handles, labels, title="Read length", loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.02), fontsize=8)
fig.suptitle("Integrated terminal read quality vs count-correlation loss, suspicious SRRs excluded", y=1.12)
fig.tight_layout(rect=(0, 0, 1, 0.95))
plt.show()


# In[13]:


srr_bins = integrated_tail[["SRR_ID", "terminal_mean_q"]].drop_duplicates().copy()
srr_bins["terminal_quality_bin"] = pd.qcut(srr_bins["terminal_mean_q"], q=4, duplicates="drop")
plot_df = integrated_tail.merge(srr_bins[["SRR_ID", "terminal_quality_bin"]], on="SRR_ID")
modes = list(pear_cols.keys())
bins = list(plot_df["terminal_quality_bin"].cat.categories)
x = np.arange(len(modes))
width = min(0.18, 0.75 / len(bins))
palette = ["#4C78A8", "#59A14F", "#F2CF5B", "#E15759"]
print(f"Data points used for plot: {len(plot_df)}")
fig, ax = plt.subplots(figsize=(12, 6.3))
for i, bin_label in enumerate(bins):
    vals = [plot_df.loc[(plot_df["terminal_quality_bin"].eq(bin_label)) & (plot_df["cleaning_mode"].eq(mode)), "correlation_loss"].dropna() for mode in modes]
    pos = x + (i - (len(bins) - 1) / 2) * width
    bp = ax.boxplot(vals, positions=pos, widths=width * 0.85, patch_artist=True, showfliers=False, manage_ticks=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(palette[i % len(palette)])
        patch.set_alpha(0.65)
    for med in bp["medians"]:
        med.set_color("#222222")
        med.set_linewidth(1.2)
    means = [v.mean() if len(v) else np.nan for v in vals]
    label = f"{bin_label.left:.1f}-{bin_label.right:.1f} (n={srr_bins['terminal_quality_bin'].eq(bin_label).sum()})"
    ax.plot(pos, means, marker="o", linewidth=1.4, color=palette[i % len(palette)], label=label)
ax.set_xticks(x)
ax.set_xticklabels(modes)
ax.set_yscale("symlog", linthresh=1e-5)
ax.set_ylim(bottom=0)
ax.set_ylabel("Correlation loss vs untrimmed counts (1 - Pearson r)")
ax.set_xlabel("Cleaning mode")
ax.set_title("Integrated terminal read quality bins vs count-correlation loss, suspicious SRRs excluded")
ax.grid(axis="y", alpha=0.25)
ax.legend(title="Terminal mean Q bin", fontsize=8)
fig.tight_layout()
plt.show()


# ### Interpretation: Integrated Quality View
# 
# The integrated quality plots ask whether samples with poorer terminal base quality are more sensitive to cleaning. After integrating MultiQC and flattened FastQC-derived quality metrics, suspicious SRRs are excluded and all remaining SRRs are plotted together rather than keeping separate MultiQC-only and FastQC-only views.
# 
# Most Adapter/P5/P10/P20 points remain close to zero correlation loss across the observed terminal quality range. Lower terminal quality does not create a clear monotonic increase in count-correlation loss for routine cleaning. The large visible deviations are concentrated in two places: P35, where stronger decay is expected, and a small non-P35 cluster near `10^0`, which is investigated separately below and treated as a library-strategy mismatch rather than ordinary RNA-seq behavior.
# 

# ### High-loss non-P35 cluster diagnostic
# 
# The integrated quality plot includes a small group of non-P35 points near a correlation loss of `10^0`. P35 is expected to behave as a stress-test condition, but Adapter/P5/P10/P20 points at this level need separate inspection. The diagnostic below identifies the affected SRRs, joins their platform metadata, and then reads the corresponding featureCounts summary files to determine whether the correlation loss reflects a meaningful cleaning effect or an unstable comparison against nearly empty count vectors.
# 

# In[14]:


import gzip
import re

high_loss_non_p35 = integrated_tail.loc[
    integrated_tail["cleaning_mode"].ne("P35") & integrated_tail["correlation_loss"].ge(0.5)
].copy()

platform_all = pd.read_csv(str(TECH / "current_srr_platform_check.tsv"), sep="\t")
high_loss_non_p35 = high_loss_non_p35.merge(
    platform_all[["SRR_ID", "project_id", "library_strategy", "library_source", "instrument_model"]],
    on=["SRR_ID", "project_id"],
    how="left",
)

high_loss_summary = (
    high_loss_non_p35
    .groupby(["project_id", "library_strategy", "library_source", "cleaning_mode"], dropna=False)
    .agg(
        rows=("SRR_ID", "size"),
        srrs=("SRR_ID", "nunique"),
        median_loss=("correlation_loss", "median"),
        min_pearson=("pearson", "min"),
        max_pearson=("pearson", "max"),
    )
    .reset_index()
    .sort_values(["project_id", "cleaning_mode"])
)

def read_featurecounts_summary(path):
    values = {}
    with gzip.open(str(path), "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "Status":
                continue
            values[parts[0]] = int(parts[1])
    total = sum(values.values())
    assigned = values.get("Assigned", 0)
    return {
        "assigned": assigned,
        "total_fragments": total,
        "assigned_pct": 100 * assigned / total if total else np.nan,
        "unmapped": values.get("Unassigned_Unmapped", 0),
        "no_features": values.get("Unassigned_NoFeatures", 0),
    }

summary_rows = []
counts_root = BASE / "flattened_counts"
mode_file = {
    "untrimmed": "untrmd_{srr}_fC.txt.summary.gz",
    "Adapter": "{srr}_trimmomatic_adapter_fC.txt.summary.gz",
    "P5": "{srr}_trimmomatic_P5_fC.txt.summary.gz",
    "P10": "{srr}_trimmomatic_P10_fC.txt.summary.gz",
    "P20": "{srr}_trimmomatic_P20_fC.txt.summary.gz",
    "P35": "{srr}_trimmomatic_P35_fC.txt.summary.gz",
}
for project_id, srr_id in high_loss_non_p35[["project_id", "SRR_ID"]].drop_duplicates().itertuples(index=False):
    for mode, pattern in mode_file.items():
        path = counts_root / project_id / pattern.format(srr=srr_id)
        if path.exists():
            row = {"project_id": project_id, "SRR_ID": srr_id, "mode": mode}
            row.update(read_featurecounts_summary(path))
            summary_rows.append(row)

high_loss_count_summary = pd.DataFrame(summary_rows)
assigned_by_mode = (
    high_loss_count_summary
    .groupby("mode", as_index=False)
    .agg(
        srrs=("SRR_ID", "nunique"),
        median_assigned=("assigned", "median"),
        median_assigned_pct=("assigned_pct", "median"),
        max_assigned_pct=("assigned_pct", "max"),
    )
    .sort_values("mode")
)

print("High-loss non-P35 rows:", len(high_loss_non_p35))
print("Affected SRRs:", high_loss_non_p35["SRR_ID"].nunique())
display(high_loss_summary)
display(assigned_by_mode)
display(
    high_loss_non_p35[
        ["SRR_ID", "project_id", "library_strategy", "cleaning_mode", "terminal_mean_q", "pearson", "correlation_loss"]
    ].sort_values("correlation_loss", ascending=False)
)


# These high-loss non-P35 points are not the manually excluded suspicious SRRs. They are all from `PRJNA1108066`, which is annotated as paired-end Illumina `miRNA-Seq`, not standard mRNA RNA-seq. The featureCounts summaries show why the correlations are unstable: the untrimmed branch assigns essentially no reads to exon-level `gene_id` features, while cleaned branches assign more reads but still remain a tiny fraction of the library. A Pearson correlation against an almost empty untrimmed vector is therefore not biologically interpretable as a routine adapter-cleaning effect.
# 
# The jump in median assigned reads for these miRNA-seq examples reflects how poorly this library type fits the exon-level mRNA counting workflow. The untrimmed branch has a median of only 2 assigned reads. Adapter/P5/P10 increase this to about 430-450 assigned reads, but the median assigned fraction is still only about 0.0014% of fragments. P20 increases the median to about 12,400 assigned reads, still only about 0.054% assigned. P35 increases the median to about 437,000 assigned reads and about 2.7% assigned, because aggressive trimming leaves much shorter fragments that can align/count somewhere in the exon annotation. In this context, shorter reads make random or nonspecific alignment more prominent rather than improving the biological quantification.
# 
# That increase should not be interpreted as improved RNA-seq quantification. It is a small-RNA library being forced through an mRNA/exon counting pipeline, with increasingly aggressive trimming making more short fragments countable. Once this miRNA-seq artifact is separated from the main interpretation, the integrated plot supports the same conclusion as the executive summary: ordinary cleaning is stable, while aggressive cleaning or mismatched library strategies can produce visible count-profile decay.
# 
# Excluding `PRJNA1108066` and excluding P35, the median correlation loss is approximately `2e-6`, the 95th percentile is about `0.0025`, the 99th percentile is about `0.0088`, and no non-P35 point remains above `0.1`. In contrast, P35 has a median correlation loss of about `0.023`, a 90th percentile near `0.89`, and 60 SRRs above `0.5`.
# 

# ## 6. Result Interpretation
# 
# The technical results support the executive interpretation: routine cleaning does not materially change gene-level count profiles for adapter-only, P5, P10, and P20 in the Illumina datasets where the mRNA-style counting workflow is appropriate. These modes remain well aligned with the untrimmed counts. If the count matrix barely changes at the technical level, downstream biological DESeq2/GSEA results are expected to be at least as stable, except in projects with unusual library behavior or workflow mismatch.
# 
# The exceptions are informative rather than contradictory. P35 shows the expected decay because it is an intentionally aggressive trimming branch. It was included as a stress-test condition, not as a routine cleaning recommendation. The non-P35 points near `10^0` are also not ordinary RNA-seq behavior; they come from `PRJNA1108066`, a miRNA-seq project that entered the broad transcriptomic filter and was forced through an mRNA/exon counting workflow. Those samples produce nearly empty untrimmed exon-level featureCounts vectors, and with shorter reads random or nonspecific alignment becomes more prominent. Their Pearson correlations should therefore be treated as a library-strategy artifact.
# 
# The biological follow-up should focus on untrimmed, adapter-only, P5, P10, and P20 in the curated RNA-seq panel. P35 should be retained only as a sensitivity/stress-test branch, and miRNA-seq or other non-standard library strategies should be excluded or flagged separately when interpreting RNA-seq cleaning sensitivity.
# 
# For endpoint downloads, the result is operational rather than universal: our recorded DB attempts do not support the older assumption that Aspera was faster and more stable in this workflow. ENA HTTP performed best in the sparse endpoint records, but the endpoint data are too unbalanced to claim that Aspera is generally inferior.
# 
