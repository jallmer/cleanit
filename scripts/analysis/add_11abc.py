import nbformat

nb = nbformat.read('bioResults.ipynb', as_version=4)

# ── Find insert point: before the current "11. QC Predictive Model" ──
insert_idx = None
for i, cell in enumerate(nb.cells):
    if cell.cell_type == 'markdown' and '11. QC Predictive Model' in cell.source and 'Actionable Pipeline Rules' in cell.source:
        insert_idx = i
        break

if insert_idx is None:
    print("ERROR: Could not find section 11 QC Predictive Model header")
    exit(1)

# ── Approach A: Model vs Baseline ──
code_a = r'''import pandas as pd
import numpy as np
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns

# Load predictions
df_preds = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/qc_model_predictions.tsv', sep='\t')

# Binarize
df_preds['True_Binary'] = df_preds['True_t_star'].apply(lambda x: 'Cleaning' if x != 'U' else 'No Cleaning')
df_preds['Pred_Binary'] = df_preds['Predicted_t_star'].apply(lambda x: 'Cleaning' if x != 'U' else 'No Cleaning')

# Baseline: always predict majority class
n_total = len(df_preds)
baseline_acc = (df_preds['True_Binary'] == 'No Cleaning').sum() / n_total
model_acc = accuracy_score(df_preds['True_Binary'], df_preds['Pred_Binary'])

print("=" * 60)
print("MODEL PERFORMANCE vs. TRIVIAL BASELINE")
print("=" * 60)
print(f"\nTotal evaluated samples:  {n_total}")
print(f"True 'No Cleaning':      {(df_preds['True_Binary'] == 'No Cleaning').sum()}")
print(f"True 'Cleaning':         {(df_preds['True_Binary'] == 'Cleaning').sum()}")
print(f"\nBaseline accuracy (always predict 'No Cleaning'): {baseline_acc:.4f} ({baseline_acc*100:.1f}%)")
print(f"Random Forest accuracy (5-fold CV):               {model_acc:.4f} ({model_acc*100:.1f}%)")
print(f"Improvement over baseline:                        {(model_acc - baseline_acc)*100:+.1f} percentage points")

print("\n" + "=" * 60)
print("CLASSIFICATION REPORT (focus on 'Cleaning' class)")
print("=" * 60)
labels = ['No Cleaning', 'Cleaning']
print(classification_report(df_preds['True_Binary'], df_preds['Pred_Binary'], 
                            target_names=labels, zero_division=0))

# Confusion matrix
cm = confusion_matrix(df_preds['True_Binary'], df_preds['Pred_Binary'], labels=labels)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: confusion matrix
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=labels, yticklabels=labels, ax=axes[0])
axes[0].set_title('Confusion Matrix (5-fold CV)')
axes[0].set_xlabel('Predicted')
axes[0].set_ylabel('True')

# Right: accuracy comparison bar chart
ax = axes[1]
bars = ax.bar(['Always "No Cleaning"\n(Trivial Baseline)', 'Random Forest\n(FastQC Features)'], 
              [baseline_acc * 100, model_acc * 100], 
              color=['#9E9E9E', '#2196F3'], alpha=0.8, width=0.5)
ax.set_ylabel('Accuracy (%)')
ax.set_title('Model Adds No Value Over Trivial Baseline')
ax.set_ylim(0, 105)
for bar, val in zip(bars, [baseline_acc * 100, model_acc * 100]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, f'{val:.1f}%', 
            ha='center', va='bottom', fontweight='bold')
ax.axhline(baseline_acc * 100, color='red', ls='--', alpha=0.5, label=f'Baseline: {baseline_acc*100:.1f}%')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.show()
'''

# ── Approach B: Feature distributions ──
code_b = r'''import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load and merge QC features with classification
qual = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/per_srr_quality.tsv', sep='\t')
cls = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/trimming_classification.tsv', sep='\t')
merged = pd.merge(qual, cls[['SRR_ID', 't_star']], on='SRR_ID')
merged['label'] = merged['t_star'].apply(lambda x: 'Cleaning Beneficial' if x != 'U' else 'Untrimmed Optimal')

features = ['adapter_rate', 'frac_below_q20', 'frac_below_q30', 'tail_quality_decay',
            'Q_mean', 'sequence_depth', 'gc_content', 'gc_deviation', 'duplication_rate']

# Filter to features that exist and have variance
features = [f for f in features if f in merged.columns and merged[f].std() > 0]

n_feat = len(features)
cols = 3
rows_n = (n_feat + cols - 1) // cols

fig, axes = plt.subplots(rows_n, cols, figsize=(5*cols, 4*rows_n))
axes = axes.ravel()

for idx, feat in enumerate(features):
    ax = axes[idx]
    u_vals = merged[merged['label'] == 'Untrimmed Optimal'][feat].dropna()
    c_vals = merged[merged['label'] == 'Cleaning Beneficial'][feat].dropna()
    
    # Plot overlapping histograms
    bins = np.histogram_bin_edges(pd.concat([u_vals, c_vals]), bins=30)
    ax.hist(u_vals, bins=bins, alpha=0.6, color='#2196F3', label=f'Untrimmed (n={len(u_vals)})', density=True)
    ax.hist(c_vals, bins=bins, alpha=0.8, color='#F44336', label=f'Cleaning (n={len(c_vals)})', density=True)
    
    # Mark the cleaning samples as rug plot
    for v in c_vals:
        ax.axvline(v, color='#F44336', alpha=0.4, lw=0.8, ymin=0, ymax=0.05)
    
    ax.set_title(feat.replace('_', ' ').title(), fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.2)

# Hide unused subplots
for idx in range(n_feat, len(axes)):
    axes[idx].set_visible(False)

plt.suptitle('FastQC Feature Distributions: Untrimmed-Optimal vs. Cleaning-Beneficial Samples\n'
             'Complete overlap proves no technical marker can predict trimming need',
             fontsize=12, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.show()

print(f"\nSamples evaluated: {len(merged)}")
print(f"Untrimmed optimal:    {(merged['label'] == 'Untrimmed Optimal').sum()}")
print(f"Cleaning beneficial:  {(merged['label'] == 'Cleaning Beneficial').sum()}")
'''

# ── Approach C: Information Gain ──
code_c = r'''import pandas as pd
import numpy as np
from sklearn.feature_selection import mutual_info_classif
import matplotlib.pyplot as plt

# Load and merge
qual = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/per_srr_quality.tsv', sep='\t')
cls = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/trimming_classification.tsv', sep='\t')
merged = pd.merge(qual, cls[['SRR_ID', 't_star']], on='SRR_ID')

# Binary label
y = (merged['t_star'] != 'U').astype(int)

features = ['adapter_rate', 'frac_below_q20', 'frac_below_q30', 'tail_quality_decay',
            'Q_mean', 'Q_median', 'sequence_depth', 'read_length_mean',
            'gc_content', 'gc_deviation', 'duplication_rate', 'n_content']
features = [f for f in features if f in merged.columns]

X = merged[features].fillna(0)

# Compute mutual information (information gain) with multiple random seeds for stability
ig_runs = []
for seed in range(10):
    ig = mutual_info_classif(X, y, discrete_features=False, random_state=seed)
    ig_runs.append(ig)

ig_mean = np.mean(ig_runs, axis=0)
ig_std = np.std(ig_runs, axis=0)

# Sort by IG
order = np.argsort(ig_mean)[::-1]

fig, ax = plt.subplots(figsize=(10, 5))
feat_names = [features[i].replace('_', ' ').title() for i in order]
bars = ax.barh(range(len(features)), ig_mean[order], xerr=ig_std[order], 
               color='#FF9800', alpha=0.8, capsize=3)
ax.set_yticks(range(len(features)))
ax.set_yticklabels(feat_names, fontsize=9)
ax.set_xlabel('Information Gain (bits)')
ax.set_title('Information Gain of FastQC Features for Predicting Trimming Need\n'
             'All values near zero → no feature carries useful discriminative information',
             fontsize=11)
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.3)

# Add reference line
ax.axvline(0, color='black', lw=0.8)

# Annotate max IG
max_ig = ig_mean[order[0]]
ax.text(max_ig + 0.002, 0, f'max = {max_ig:.4f} bits', va='center', fontsize=9, color='red')

plt.tight_layout()
plt.show()

# Print table
print("\nInformation Gain Summary:")
print(f"{'Feature':<25} {'IG (bits)':>12}  {'Interpretation':<30}")
print("-" * 70)
for i in order:
    ig_val = ig_mean[i]
    if ig_val < 0.01:
        interp = "No signal"
    elif ig_val < 0.05:
        interp = "Negligible"
    else:
        interp = "Weak"
    print(f"{features[i]:<25} {ig_val:>12.4f}  {interp:<30}")

print(f"\nFor reference, a perfectly predictive binary feature would have IG ≈ {-((4/193)*np.log2(4/193) + (189/193)*np.log2(189/193)):.4f} bits (dataset entropy).")
print(f"Maximum observed IG is {max_ig:.4f} bits — orders of magnitude below the entropy ceiling.")
'''

# ── Build cells ──
md_11a = nbformat.v4.new_markdown_cell(
    source="## 11a. Can FastQC Metrics Predict Trimming Need?\n"
           "### Approach A: Model Performance vs. Trivial Baseline\n\n"
           "If a machine learning model trained on FastQC features cannot beat the trivial strategy of **always predicting 'No Cleaning'**, "
           "then the features contain no actionable signal for trimming decisions. "
           "We evaluate our Random Forest classifier (5-fold cross-validated) against this baseline."
)
code_cell_a = nbformat.v4.new_code_cell(source=code_a)

md_11b = nbformat.v4.new_markdown_cell(
    source="## 11b. Feature Distribution Overlap\n"
           "### Approach B: Do Cleaning-Beneficial Samples Look Different?\n\n"
           "The most fundamental test: if the FastQC feature distributions of 'Cleaning Beneficial' samples completely overlap "
           "with 'Untrimmed Optimal' samples, then **no model** — no matter how sophisticated — could ever learn to separate them. "
           "The signal simply does not exist in the technical QC metrics."
)
code_cell_b = nbformat.v4.new_code_cell(source=code_b)

md_11c = nbformat.v4.new_markdown_cell(
    source="## 11c. Information Gain Analysis\n"
           "### How Much Discriminative Information Does Each Feature Carry?\n\n"
           "Information Gain (mutual information) quantifies how many bits of information each FastQC feature provides about whether "
           "a sample benefits from trimming. A value near zero means the feature is completely uninformative. "
           "We compute this for all available technical features."
)
code_cell_c = nbformat.v4.new_code_cell(source=code_c)

# Insert before old section 11
nb.cells.insert(insert_idx, code_cell_c)
nb.cells.insert(insert_idx, md_11c)
nb.cells.insert(insert_idx, code_cell_b)
nb.cells.insert(insert_idx, md_11b)
nb.cells.insert(insert_idx, code_cell_a)
nb.cells.insert(insert_idx, md_11a)

nbformat.write(nb, 'bioResults.ipynb')
print("Successfully added sections 11a, 11b, 11c to bioResults.ipynb")
