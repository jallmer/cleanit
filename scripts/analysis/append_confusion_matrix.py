import nbformat

nb_path = "bioResults.ipynb"
nb = nbformat.read(nb_path, as_version=4)

new_markdown = """## 11. QC Prediction Model (Technical FastQC Features vs Biological Need)
**Interpretation:** This section evaluates whether standard FastQC technical metrics (e.g. Q-scores, sequence depth, GC content) can reliably predict if a sample needs trimming based on our strict biological Jaccard criteria.

We trained a Random Forest classifier using 5-fold cross-validation.
*   **No Cleaning:** The optimal biological choice was to leave the data untrimmed (`U`).
*   **Cleaning Needed:** The optimal biological choice was a trimming method (e.g., `P20`, `P35`).

A strong diagonal indicates that technical metrics perfectly predict biological necessity. Off-diagonal elements indicate where FastQC metrics might suggest trimming when it's biologically harmful, or vice versa."""

new_code = """import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

df_preds = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/qc_model_predictions.tsv', sep='\\t')

# Binarize the classes
# True: actual empirical biological need
df_preds['True_Binary'] = df_preds['True_t_star'].apply(lambda x: 'No Cleaning' if x == 'U' else 'Cleaning Needed')

# Predicted: what the Random Forest predicted
df_preds['Predicted_Binary'] = df_preds['Predicted_t_star'].apply(lambda x: 'No Cleaning' if x == 'U' else 'Cleaning Needed')

labels = ['No Cleaning', 'Cleaning Needed']
cm = confusion_matrix(df_preds['True_Binary'], df_preds['Predicted_Binary'], labels=labels)

plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=labels, yticklabels=labels)
plt.title('Confusion Matrix: FastQC Metrics vs Biological Need')
plt.xlabel('Predicted by Random Forest (FastQC Features)')
plt.ylabel('True Optimal Biological Choice (t*)')
plt.tight_layout()
plt.show()"""

md_cell = nbformat.v4.new_markdown_cell(new_markdown)
code_cell = nbformat.v4.new_code_cell(new_code)

# Avoid adding it multiple times
has_cm = any('Confusion Matrix: FastQC' in cell.source for cell in nb.cells)
if not has_cm:
    nb.cells.extend([md_cell, code_cell])
    nbformat.write(nb, nb_path)
    print("Appended Confusion Matrix to bioResults.ipynb")
else:
    print("Confusion Matrix already in bioResults.ipynb")
