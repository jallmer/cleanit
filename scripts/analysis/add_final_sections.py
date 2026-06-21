import nbformat

def add_final_sections(nb_path):
    with open(nb_path) as f:
        nb = nbformat.read(f, as_version=4)

    # Remove existing final sections if they exist to allow easy re-running
    cells_to_keep = []
    for cell in nb.cells:
        if "10. Gene vs Pathway Robustness" in cell.source or "11. QC Predictive Model" in cell.source:
            break
        cells_to_keep.append(cell)
    nb.cells = cells_to_keep

    md_10 = nbformat.v4.new_markdown_cell("""## 10. Gene vs Pathway Robustness
### Interpretation & Biological Significance

This section provides the fundamental justification for using **Gene Set Enrichment Analysis (GSEA)** rather than relying solely on raw lists of Differentially Expressed Genes (DEGs) to evaluate bioinformatics pipeline stability.

*   **The Extreme Volatility of DEGs (Red Boxes):** Notice how the Jaccard similarity for DEGs (the overlap of strictly significant genes) is almost entirely squashed below `0.1`! This indicates that simply dropping a single sample from the analysis (Leave-One-Out) completely destroys the stability of the "significant" gene list. The resulting lists of genes have almost zero overlap. Relying on strict p-value cutoffs for single genes is incredibly sensitive to statistical noise and sample variance.
*   **The Stability of Pathways (Green Boxes):** Conversely, the GSEA Pathway Jaccard similarity remains tightly clustered above `0.9` for all reasonable trimming methods. Even when the specific underlying genes fluctuate wildly, the overarching biological pathways they map to remain incredibly stable. (Note: `P35` is the exception, where the extremely aggressive trimming destroys so much data that even pathway concordance begins to fail).

**Conclusion:** When evaluating whether a trimming method "damages" or "benefits" an experiment, looking at individual genes provides too much noise. Pathway-level concordance proves that the core biological narrative of the experiment is robust, making it the superior metric for our trimming classification.""")

    code_10 = nbformat.v4.new_code_cell("""import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load all concordance data
df_all = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/all_concordance.tsv", sep="\\t")

methods = ["U", "A", "P5", "P10", "P20", "P35"]

# Prepare data for matplotlib side-by-side boxplot
gene_data = [df_all[f"{m}_jaccard_deg"].dropna().values for m in methods]
path_data = [df_all[f"{m}_jaccard_path"].dropna().values for m in methods]

fig, ax = plt.subplots(figsize=(12, 6))

# Plot Gene Level
bplot1 = ax.boxplot(gene_data, positions=np.arange(len(methods)) * 2.0 - 0.4, widths=0.6, patch_artist=True, tick_labels=["" for _ in methods])
for patch in bplot1['boxes']:
    patch.set_facecolor('#F44336')

# Plot Pathway Level
bplot2 = ax.boxplot(path_data, positions=np.arange(len(methods)) * 2.0 + 0.4, widths=0.6, patch_artist=True, tick_labels=["" for _ in methods])
for patch in bplot2['boxes']:
    patch.set_facecolor('#4CAF50')

# Custom legend and labels
ax.plot([], [], color='#F44336', marker='s', markersize=10, linestyle='None', label='Strict DEG List (p < 0.05)')
ax.plot([], [], color='#4CAF50', marker='s', markersize=10, linestyle='None', label='GSEA Pathways (p < 0.05)')

ax.set_xticks(np.arange(len(methods)) * 2.0)
ax.set_xticklabels(methods)

plt.title("Concordance Stability: Gene vs Pathway Level")
plt.ylabel("Jaccard Similarity to Consensus Reference")
plt.ylim(0, 1.05)
plt.grid(axis='y', alpha=0.3)
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()
""")

    md_11 = nbformat.v4.new_markdown_cell("""## 11. QC Predictive Model
### Interpretation & Actionable Pipeline Rules

Having proven that pathway-level analysis is robust, we utilized it to score the biological "benefit" of different trimming methods ($t^*$) across **809 fully evaluated samples**. We mapped these optimal trimming choices back to the raw, pre-alignment **FastQC** metrics.

1.  **Feature Importance (Left Plot):** This plot reveals which raw sequencing quality metrics are the most predictive of trimming benefit. Metrics like `frac_below_q20` (the percentage of low-quality bases), `sequence_depth`, and `adapter_rate` dictate whether trimming will help or harm the downstream biological conclusions.
2.  **Model Coefficients (Right Plot):** This plot visualizes the directional impact of each metric.
    *   **Positive Coefficients** (e.g., high adapter contamination, high low-quality base fractions) strongly push the model to recommend aggressive trimming (like `A` or `P10`).
    *   **Negative Coefficients** (e.g., standard high depth, normal GC content) push the model to recommend **Untrimmed (U)**. As we observed in the classification stage, if trimming does not strictly improve the pathway concordance, the biologically superior choice is to leave the data untrimmed to avoid introducing alignment bias.

**Final Pipeline Recommendation:** Because "Neutral" trimming methods provide no biological benefit but risk mathematical bias, our classification heavily favors "Untrimmed" for the vast majority of modern high-quality RNA-seq runs. Trimming is dynamically recommended *only* when the predictive FastQC metrics (like extreme adapter bleed-through) trigger the model's threshold.""")

    code_11 = nbformat.v4.new_code_cell("""import pandas as pd
from IPython.display import Image, display

try:
    print("QC Model Predictors:")
    display(Image(filename="/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/plots/qc_predictors.png"))
    
    print("\\nQC Model Coefficients:")
    display(Image(filename="/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/plots/qc_coefficients.png"))
    
except Exception as e:
    print(f"Error loading QC images: {e}")
""")

    md_10_1 = nbformat.v4.new_markdown_cell("""## 10.1 Gene Rank Correlation Heatmap
**Interpretation:** This heatmap visualizes the Spearman correlation (`rho_gene`) between the Candidate DEGs and the Reference DEGs. Aggressive methods like `P35` consistently destroy gene rank correlation.""")

    code_10_1 = nbformat.v4.new_code_cell("""import seaborn as sns
# Pivot data for heatmap
methods_order = ["U", "A", "P5", "P10", "P20", "P35"]
heat_data = pd.DataFrame(index=df_all['project'].unique(), columns=methods_order)
for p in df_all['project'].unique():
    subset = df_all[df_all['project'] == p]
    for m in methods_order:
        heat_data.loc[p, m] = subset[f"{m}_rho_gene"].mean()

heat_data = heat_data.apply(pd.to_numeric)
plt.figure(figsize=(10, 10))
sns.heatmap(heat_data, annot=True, cmap="viridis", fmt=".2f")
plt.title("Mean Spearman Correlation (rho_gene) per Project")
plt.ylabel("Project")
plt.xlabel("Trimming Method")
plt.show()""")

    md_10_2 = nbformat.v4.new_markdown_cell("""## 10.2 DEG Count Scatter Plot
**Interpretation:** This scatter plot compares the number of Significant DEGs in the Reference vs the Candidate. Notice how `P35` (brown/purple dots) skyrockets above the $y=x$ line, indicating massive false-positive inflation.""")

    code_10_2 = nbformat.v4.new_code_cell("""fig, ax = plt.subplots(figsize=(8, 8))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

for i, m in enumerate(methods_order):
    ref_col = f"{m}_n_ref_degs"
    cand_col = f"{m}_n_cand_degs"
    ax.scatter(df_all[ref_col], df_all[cand_col], label=m, alpha=0.5, s=20, c=colors[i])

# y=x reference line
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  
    np.max([ax.get_xlim(), ax.get_ylim()]),  
]
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_xlim(lims)
ax.set_ylim(lims)

ax.set_xlabel("Number of DEGs (Reference)")
ax.set_ylabel("Number of DEGs (Candidate)")
ax.set_title("DEG Count: Reference vs Candidate")
ax.legend()
plt.show()""")

    md_10_3 = nbformat.v4.new_markdown_cell("""## 10.3 Pathway Concordance Dashboard
**Interpretation:** This dashboard provides a deeper dive into the pathway-level metrics. It shows the Jaccard overlap, the Spearman rank correlation of Normalized Enrichment Scores (NES), and the Direction Concordance of the top 50 pathways.""")

    code_10_3 = nbformat.v4.new_code_cell("""fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# 1. Jaccard Pathway
path_jaccard_data = [df_all[f"{m}_jaccard_path"].dropna().values for m in methods_order]
axes[0].boxplot(path_jaccard_data, tick_labels=methods_order, patch_artist=True)
axes[0].set_title("Pathway Jaccard Overlap")
axes[0].set_ylim(0, 1.05)

# 2. Rho Pathway
path_rho_data = [df_all[f"{m}_rho_path"].dropna().values for m in methods_order]
axes[1].boxplot(path_rho_data, tick_labels=methods_order, patch_artist=True)
axes[1].set_title("Pathway Spearman Correlation (rho)")
axes[1].set_ylim(-1.05, 1.05)

# 3. Direction Concordance
path_dir_data = [df_all[f"{m}_dir_concordance"].dropna().values for m in methods_order]
axes[2].boxplot(path_dir_data, tick_labels=methods_order, patch_artist=True)
axes[2].set_title("Top 50 Pathway Direction Concordance")
axes[2].set_ylim(0, 1.05)

for ax in axes:
    for patch in ax.patches:
        patch.set_facecolor('#4CAF50')

plt.tight_layout()
plt.show()""")

    nb.cells.extend([md_10, code_10, md_10_1, code_10_1, md_10_2, code_10_2, md_10_3, code_10_3, md_11, code_11])
    with open(nb_path, "w") as f:
        nbformat.write(nb, f)
    print("Successfully added final sections to notebook.")

if __name__ == "__main__":
    add_final_sections("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/bioResults.ipynb")
