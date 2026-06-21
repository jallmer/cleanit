import nbformat

def add_cells(nb_path):
    with open(nb_path) as f:
        nb = nbformat.read(f, as_version=4)

    # Check if section 9 already exists
    if any("9. Sample-Specific Trimming Classification" in cell.source for cell in nb.cells):
        print("Cells already exist.")
        return

    md_cell = nbformat.v4.new_markdown_cell("""## 9. Sample-Specific Trimming Classification\nThe previous sections evaluated trimming on a whole-project basis. This section uses the rigorous Leave-One-Sample-Out (LOSO) methodology to classify whether trimming a *specific* sample is Helpful, Neutral, or Harmful to the project's biological conclusion, after penalizing for technical information loss.""")

    code_cell = nbformat.v4.new_code_cell("""import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load classification results
cls_file = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/trimming_classification.tsv"

try:
    df_cls = pd.read_csv(cls_file, sep="\\t")
    
    # Count classifications per method
    methods = ["A", "P5", "P10", "P20", "P35"]
    counts = {"helpful": [], "neutral": [], "harmful": []}
    
    for m in methods:
        c = df_cls[f"{m}_class"].value_counts()
        for status in counts.keys():
            counts[status].append(c.get(status, 0))
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Classification plot
    ax = axes[0]
    bottom = np.zeros(len(methods))
    colors = {"helpful": "#4CAF50", "neutral": "#9E9E9E", "harmful": "#F44336"}
    
    for status, values in counts.items():
        ax.bar(methods, values, label=status, bottom=bottom, color=colors[status], alpha=0.8)
        bottom += np.array(values)
    
    ax.set_title("Sample-Specific Trimming Classification")
    ax.set_ylabel("Number of Samples")
    ax.set_xlabel("Trimming Method")
    ax.legend(title="Classification")
    ax.grid(axis="y", alpha=0.3)
    
    # Optimal method plot
    ax2 = axes[1]
    t_star_counts = df_cls["t_star"].value_counts().reindex(["U"] + methods).fillna(0)
    ax2.bar(t_star_counts.index, t_star_counts.values, color="#2196F3", alpha=0.8)
    ax2.set_title("Optimal Trimming Choice (t*) per Sample")
    ax2.set_ylabel("Number of Samples")
    ax2.set_xlabel("Preferred Method")
    ax2.grid(axis="y", alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    print(f"Evaluated {len(df_cls)} sample(s) using Leave-One-Sample-Out.")
except FileNotFoundError:
    print("Run 06_run_de_gsea.py and 07_classify_trimming.py first to generate the classification dataset.")
""")

    nb.cells.extend([md_cell, code_cell])
    with open(nb_path, "w") as f:
        nbformat.write(nb, f)
    print("Successfully added new section to notebook.")

add_cells("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/bioResults.ipynb")
