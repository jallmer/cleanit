#!/usr/bin/env python3
"""
Retrain the QC prediction model with the full per_srr_quality.tsv (830 SRRs).
Generates qc_model_predictions.tsv for use in bioResults.ipynb §11a.
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import LabelEncoder

ANALYSIS = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis"

# Load data
qual = pd.read_csv(f"{ANALYSIS}/per_srr_quality.tsv", sep="\t")
cls = pd.read_csv(f"{ANALYSIS}/trimming_classification.tsv", sep="\t")

# Merge
merged = pd.merge(cls[["SRR_ID", "project_id", "t_star"]], qual, on="SRR_ID")
print(f"Merged samples: {len(merged)}")
print(f"t_star distribution: {merged['t_star'].value_counts().to_dict()}")

# Features (use all available, drop NaN rows per feature set)
features_full = ["Q_mean", "Q_median", "sequence_depth", "read_length_mean",
                 "tail_quality_decay", "gc_content", "gc_deviation",
                 "adapter_rate", "frac_below_q20", "frac_below_q30",
                 "duplication_rate", "n_content"]

# Drop rows with any NaN in features
available = merged.dropna(subset=features_full).copy()
print(f"Samples with all features: {len(available)}")

X = available[features_full].values
y = available["t_star"].values

# 5-fold stratified CV
# Since we have very few positive cases, use stratified folds
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
clf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight="balanced")

y_pred = cross_val_predict(clf, X, y, cv=cv)

# Save predictions
preds = pd.DataFrame({
    "SRR_ID": available["SRR_ID"].values,
    "project_id": available["project_id_x"].values if "project_id_x" in available.columns else available["project_id"].values,
    "True_t_star": y,
    "Predicted_t_star": y_pred,
})
preds.to_csv(f"{ANALYSIS}/concordance/qc_model_predictions.tsv", sep="\t", index=False)
print(f"\nSaved predictions for {len(preds)} samples to qc_model_predictions.tsv")

# Quick summary
from sklearn.metrics import accuracy_score, classification_report
print(f"\nAccuracy: {accuracy_score(y, y_pred):.4f}")
print(f"Baseline (always U): {(y == 'U').sum() / len(y):.4f}")
print(f"\nPrediction distribution: {pd.Series(y_pred).value_counts().to_dict()}")
print(f"True distribution: {pd.Series(y).value_counts().to_dict()}")
