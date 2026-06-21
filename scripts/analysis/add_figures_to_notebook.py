import nbformat

nb_path = "bioResults.ipynb"
nb = nbformat.read(nb_path, as_version=4)

markdown_cell = nbformat.v4.new_markdown_cell("""## Whole-Project Concordance

We now evaluate the whole-project concordance (comparing the full trimmed dataset against the untrimmed reference).
The concordance table contains 40 projects and reports DEG Jaccard index, log2FC correlation, and GSEA NES correlation.
""")

code_cell = nbformat.v4.new_code_cell("""import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

CONCORDANCE_DIR = Path("concordance")
df = pd.read_csv(CONCORDANCE_DIR / "whole_project_concordance.tsv", sep="\\t")
METHODS = ["A", "P5", "P10", "P20", "P35"]

# Set up styling
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
palette = sns.color_palette("husl", len(METHODS))

# 1. DEG Jaccard Overlap
plt.figure(figsize=(8, 5))
sns.boxplot(data=df, x="method", y="deg_jaccard", order=METHODS, palette=palette)
sns.stripplot(data=df, x="method", y="deg_jaccard", order=METHODS, color=".3", size=4, alpha=0.6)
plt.title("DEG Jaccard Overlap (Trimmed vs Untrimmed)")
plt.ylabel("Jaccard Index")
plt.xlabel("Trimming Method")
plt.ylim(0, 1.05)
plt.show()

# 2. DEG log2FC Correlation
plt.figure(figsize=(8, 5))
sns.boxplot(data=df, x="method", y="rho_log2fc", order=METHODS, palette=palette)
sns.stripplot(data=df, x="method", y="rho_log2fc", order=METHODS, color=".3", size=4, alpha=0.6)
plt.title("Spearman Correlation of log2FC (Shared DEGs)")
plt.ylabel(r"Spearman $\\rho$")
plt.xlabel("Trimming Method")
plt.ylim(0, 1.05)
plt.show()

# 3. GSEA NES Correlation
plt.figure(figsize=(8, 5))
sns.boxplot(data=df, x="method", y="rho_nes", order=METHODS, palette=palette)
sns.stripplot(data=df, x="method", y="rho_nes", order=METHODS, color=".3", size=4, alpha=0.6)
plt.title("GSEA NES Correlation (Shared Pathways)")
plt.ylabel(r"Spearman $\\rho$")
plt.xlabel("Trimming Method")
plt.ylim(0, 1.05)
plt.show()

# 4. GSEA Directional Concordance (Top 20)
plt.figure(figsize=(8, 5))
sns.boxplot(data=df, x="method", y="dir_concordance", order=METHODS, palette=palette)
sns.stripplot(data=df, x="method", y="dir_concordance", order=METHODS, color=".3", size=4, alpha=0.6)
plt.title("GSEA Directional Concordance (Top 20 Pathways)")
plt.ylabel("Proportion of Same Sign NES")
plt.xlabel("Trimming Method")
plt.ylim(0, 1.05)
plt.show()
""")

nb.cells.extend([markdown_cell, code_cell])
nbformat.write(nb, nb_path)
print(f"Appended figure generation cells to {nb_path}")
