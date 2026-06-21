import nbformat

nb = nbformat.read('bioResults.ipynb', as_version=4)

# ── Cell code for 9c: Outlier Diagnostic Table ──
outlier_code = r'''import pandas as pd
from IPython.display import display, HTML

cls = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/trimming_classification.tsv", sep="\t")
outliers = cls[cls["t_star"] != "U"].copy()

bio = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/bio_concordance.tsv", sep="\t")
piv_deg = bio.pivot_table(index="SRR_ID", columns="method", values="jaccard_deg")
piv_path = bio.pivot_table(index="SRR_ID", columns="method", values="jaccard_pathway")

trim = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/trimmomatic_detail.tsv", sep="\t")
mode_map = {"A": "adapter_only", "P5": "P5", "P10": "P10", "P20": "P20", "P35": "P35"}

rows = []
for _, row in outliers.iterrows():
    srr, proj, tstar = row["SRR_ID"], row["project_id"], row["t_star"]
    deg_u  = piv_deg.loc[srr, "U"]    if srr in piv_deg.index  else None
    deg_t  = piv_deg.loc[srr, tstar]  if srr in piv_deg.index  else None
    path_u = piv_path.loc[srr, "U"]   if srr in piv_path.index else None
    path_t = piv_path.loc[srr, tstar] if srr in piv_path.index else None
    delta_deg  = deg_t  - deg_u  if None not in (deg_u, deg_t)   else None
    delta_path = path_t - path_u if None not in (path_u, path_t) else None

    trim_sub = trim[(trim["SRR_ID"] == srr) & (trim["mode"] == mode_map.get(tstar, tstar))]
    lost = trim_sub["dropped_pct"].values[0] if len(trim_sub) > 0 else None

    rows.append({
        "SRR": srr, "Project": proj, "t*": tstar,
        "Δ JD_gene": round(delta_deg, 4) if delta_deg is not None else None,
        "Δ JD_path": round(delta_path, 4) if delta_path is not None else None,
        "Reads Lost (%)": round(lost, 2) if lost is not None else None,
    })

df_out = pd.DataFrame(rows)
# Highlight project clusters
display(df_out.style
    .format(precision=4, na_rep="—")
    .set_caption("Outlier samples where t* ≠ Untrimmed (strict gene+pathway concordance)")
    .set_table_styles([{"selector": "caption", "props": [("font-weight", "bold"), ("font-size", "13px")]}])
)

# Summary
n_total = len(cls)
n_out = len(outliers)
projects = outliers["project_id"].value_counts()
print(f"\n{n_out}/{n_total} samples ({100*n_out/n_total:.1f}%) prefer trimming.")
print(f"These come from {len(projects)} projects.  Project breakdown:")
for p, c in projects.items():
    print(f"  {p}: {c} sample(s)")
print(f"\nMedian Δ JD_gene:  {df_out['Δ JD_gene'].median():.4f}")
print(f"Median Δ JD_path:  {df_out['Δ JD_path'].median():.4f}")
print(f"Median reads lost: {df_out['Reads Lost (%)'].median():.2f}%")
'''

# ── Cell code for 9d: Read-loss vs Benefit scatter ──
scatter_code = r'''import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

cls = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/trimming_classification.tsv", sep="\t")
bio = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/bio_concordance.tsv", sep="\t")
trim = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/trimmomatic_detail.tsv", sep="\t")

piv_deg = bio.pivot_table(index="SRR_ID", columns="method", values="jaccard_deg")
piv_path = bio.pivot_table(index="SRR_ID", columns="method", values="jaccard_pathway")
mode_map = {"A": "adapter_only", "P5": "P5", "P10": "P10", "P20": "P20", "P35": "P35"}
methods = ["A", "P5", "P10", "P20", "P35"]

rows = []
for srr in cls["SRR_ID"]:
    if srr not in piv_deg.index:
        continue
    u_deg = piv_deg.loc[srr, "U"]
    u_path = piv_path.loc[srr, "U"] if srr in piv_path.index else None

    # For each method, record its delta and read loss
    for m in methods:
        if m not in piv_deg.columns or pd.isna(piv_deg.loc[srr, m]):
            continue
        delta_deg = piv_deg.loc[srr, m] - u_deg
        delta_path = (piv_path.loc[srr, m] - u_path) if u_path is not None and m in piv_path.columns and pd.notnull(piv_path.loc[srr, m]) else None
        
        trim_sub = trim[(trim["SRR_ID"] == srr) & (trim["mode"] == mode_map.get(m, m))]
        lost = trim_sub["dropped_pct"].values[0] if len(trim_sub) > 0 else None
        
        if lost is not None:
            rows.append({
                "SRR_ID": srr,
                "method": m,
                "delta_deg": delta_deg,
                "delta_path": delta_path,
                "reads_lost_pct": lost,
            })

df = pd.DataFrame(rows)

# Color by method
colors = {"A": "#4CAF50", "P5": "#2196F3", "P10": "#FF9800", "P20": "#9C27B0", "P35": "#F44336"}

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left: Gene Jaccard delta vs reads lost
ax = axes[0]
for m in methods:
    sub = df[df["method"] == m]
    ax.scatter(sub["reads_lost_pct"], sub["delta_deg"], c=colors[m], alpha=0.4, s=15, label=m, edgecolors="none")
ax.axhline(0.01, color="green", ls="--", alpha=0.6, label="Threshold (+0.01)")
ax.axhline(-0.01, color="red", ls="--", alpha=0.6, label="Threshold (−0.01)")
ax.axhline(0, color="black", ls="-", alpha=0.3)
ax.set_xlabel("Reads Lost (%)")
ax.set_ylabel("Δ Gene Jaccard (method − Untrimmed)")
ax.set_title("Gene Jaccard Improvement vs. Read Cost")
ax.legend(fontsize=8, loc="lower left")
ax.grid(alpha=0.2)

# Right: Pathway Jaccard delta vs reads lost
ax = axes[1]
for m in methods:
    sub = df[df["method"] == m].dropna(subset=["delta_path"])
    ax.scatter(sub["reads_lost_pct"], sub["delta_path"], c=colors[m], alpha=0.4, s=15, label=m, edgecolors="none")
ax.axhline(0.01, color="green", ls="--", alpha=0.6, label="Threshold (+0.01)")
ax.axhline(-0.01, color="red", ls="--", alpha=0.6, label="Threshold (−0.01)")
ax.axhline(0, color="black", ls="-", alpha=0.3)
ax.set_xlabel("Reads Lost (%)")
ax.set_ylabel("Δ Pathway Jaccard (method − Untrimmed)")
ax.set_title("Pathway Jaccard Improvement vs. Read Cost")
ax.legend(fontsize=8, loc="lower left")
ax.grid(alpha=0.2)

plt.suptitle("Trade-off: Biological Benefit vs. Information Loss from Trimming", fontsize=13, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# Summary stats
print("For mild methods (A/P5/P10):")
mild = df[df["method"].isin(["A", "P5", "P10"])]
print(f"  Reads lost: median {mild['reads_lost_pct'].median():.2f}%, max {mild['reads_lost_pct'].max():.2f}%")
print(f"  Δ Gene JD:  median {mild['delta_deg'].median():.4f}")
print(f"  Δ Path JD:  median {mild['delta_path'].dropna().median():.4f}")
print(f"\nFor P20:")
p20 = df[df["method"] == "P20"]
print(f"  Reads lost: median {p20['reads_lost_pct'].median():.2f}%, max {p20['reads_lost_pct'].max():.2f}%")
print(f"  Δ Gene JD:  median {p20['delta_deg'].median():.4f}")
print(f"  Δ Path JD:  median {p20['delta_path'].dropna().median():.4f}")
print(f"\nFor P35:")
p35 = df[df["method"] == "P35"]
print(f"  Reads lost: median {p35['reads_lost_pct'].median():.2f}%, max {p35['reads_lost_pct'].max():.2f}%")
print(f"  Δ Gene JD:  median {p35['delta_deg'].median():.4f}")
print(f"  Δ Path JD:  median {p35['delta_path'].dropna().median():.4f}")
'''

# ── Find where to insert: after section 9b ──
insert_idx = None
for i, cell in enumerate(nb.cells):
    if cell.cell_type == 'markdown' and '10. Gene vs Pathway Robustness' in cell.source:
        insert_idx = i
        break

if insert_idx is None:
    print("ERROR: Could not find section 10 header")
    exit(1)

# Build the new cells
md_9c = nbformat.v4.new_markdown_cell(
    source="## 9c. Outlier Diagnostic: Samples Where Trimming Appears Beneficial\n"
           "Of the 830 samples evaluated, only a small minority have `t* ≠ U` under the strict concordance criterion. "
           "This table examines those outliers to determine whether the benefit is biologically meaningful or an artefact of marginal threshold-crossing."
)
code_9c = nbformat.v4.new_code_cell(source=outlier_code)

md_9d = nbformat.v4.new_markdown_cell(
    source="## 9d. Trade-off: Biological Benefit vs. Information Loss\n"
           "For every sample × method combination, we plot the Jaccard improvement (Δ) against the percentage of reads lost to trimming. "
           "If trimming were genuinely beneficial, we would expect points in the upper-left quadrant (large improvement, low read loss). "
           "Instead, the vast majority of points cluster around Δ ≈ 0, while aggressive methods (P20, P35) discard substantial fractions of reads for no measurable biological gain."
)
code_9d = nbformat.v4.new_code_cell(source=scatter_code)

# Insert before section 10
nb.cells.insert(insert_idx, code_9d)
nb.cells.insert(insert_idx, md_9d)
nb.cells.insert(insert_idx, code_9c)
nb.cells.insert(insert_idx, md_9c)

nbformat.write(nb, 'bioResults.ipynb')
print("Successfully added sections 9c and 9d to bioResults.ipynb")
