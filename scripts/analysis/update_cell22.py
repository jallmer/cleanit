import nbformat

nb_path = "bioResults.ipynb"
nb = nbformat.read(nb_path, as_version=4)

new_source = """import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load LONG-format concordance data for per-sample plots
df_all = pd.read_csv("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/bio_concordance.tsv", sep="\\t")

methods = ["U", "A", "P5", "P10", "P20", "P35"]

# Prepare data for matplotlib side-by-side boxplot
gene_data = [df_all[df_all["method"] == m]["jaccard_deg"].dropna().values for m in methods]
path_data = [df_all[df_all["method"] == m]["jaccard_pathway"].dropna().values for m in methods]

fig, ax = plt.subplots(figsize=(12, 6))

# Plot Gene Level
bplot1 = ax.boxplot(gene_data, positions=np.arange(len(methods)) * 2.0 - 0.4, widths=0.6, patch_artist=True, tick_labels=["" for _ in methods])
for patch in bplot1["boxes"]:
    patch.set_facecolor("#F44336")

# Plot Pathway Level
bplot2 = ax.boxplot(path_data, positions=np.arange(len(methods)) * 2.0 + 0.4, widths=0.6, patch_artist=True, tick_labels=["" for _ in methods])
for patch in bplot2["boxes"]:
    patch.set_facecolor("#4CAF50")

# Custom legend and labels
ax.plot([], [], color="#F44336", marker="s", markersize=10, linestyle="None", label="Strict DEG List (p < 0.05)")
ax.plot([], [], color="#4CAF50", marker="s", markersize=10, linestyle="None", label="GSEA Pathways (p < 0.05)")

ax.set_xticks(np.arange(len(methods)) * 2.0)
ax.set_xticklabels(methods)

plt.title("Concordance Stability: Gene vs Pathway Level")
plt.ylabel("Jaccard Similarity to Consensus Reference")
plt.ylim(0, 1.05)
plt.grid(axis="y", alpha=0.3)
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()"""

# Find Cell 22 and update it
for cell in nb.cells:
    if cell.cell_type == "code" and "all_concordance.tsv" in cell.source:
        cell.source = new_source
        print("Updated Cell 22 to use bio_concordance.tsv.")
        break

nbformat.write(nb, nb_path)
print("Saved bioResults.ipynb.")
