#!/usr/bin/env python3
"""
08_fit_qc_model.py — model QC predictors of trimming benefit and aggressiveness.

Uses the merged sample feature table so the model can consume:
  - raw FastQC predictors
  - t* / benefit scores from the concordance step
  - technical penalty features for the chosen and candidate methods
"""

from __future__ import annotations

import csv
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")

OUT_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis"
PLOTS_DIR = os.path.join(OUT_DIR, "plots")
FEATURE_TABLE = os.path.join(OUT_DIR, "sample_feature_table.tsv")

METHOD_ORDINAL = {"U": 0, "A": 1, "P5": 2, "P10": 3, "P20": 4, "P35": 5}
QC_FEATURES = [
    "read_length_mean",
    "sequence_depth",
    "Q_mean",
    "frac_below_q20",
    "frac_below_q30",
    "tail_quality_decay",
    "adapter_rate",
    "duplication_rate",
    "n_content",
    "gc_deviation",
]


def safe_float(v):
    if v is None or v == "" or str(v).strip() == "NA":
        return np.nan
    try:
        return float(str(v).rstrip("%"))
    except ValueError:
        return np.nan


def load_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_simple_table(path, rows):
    if not rows:
        return
    cols = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in cols) + "\n")


def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)
    rows = load_tsv(FEATURE_TABLE)
    if not rows:
        print(f"ERROR: {FEATURE_TABLE} not found or empty. Run 13_build_sample_feature_table.py first.")
        return

    merged = []
    for row in rows:
        benefit = safe_float(row.get("benefit_B"))
        benefit_net = safe_float(row.get("benefit_B_net"))
        t_star = row.get("t_star", "U")
        entry = {
            "SRR_ID": row.get("SRR_ID", ""),
            "project_id": row.get("project_id", ""),
            "benefit_B": benefit,
            "benefit_B_net": benefit_net,
            "t_star": t_star,
            "t_star_ordinal": METHOD_ORDINAL.get(t_star, 0),
        }
        for feature in QC_FEATURES:
            entry[feature] = safe_float(row.get(feature))
        merged.append(entry)

    print(f"Loaded {len(merged)} SRRs from sample feature table")

    try:
        import pandas as pd
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
        from statsmodels.miscmodels.ordinal_model import OrderedModel
    except ImportError:
        print("WARNING: pandas/statsmodels not available. Skipping model fitting.")
        return

    df = pd.DataFrame(merged)
    df["project_id"] = df["project_id"].astype(str)

    coeff_rows = []
    clean = df.dropna(subset=QC_FEATURES + ["benefit_B"])
    print(f"Rows with complete QC + benefit_B data: {len(clean)}")
    if len(clean) >= 10:
        for feature in QC_FEATURES:
            std = clean[feature].std()
            if std and std > 0:
                clean.loc[:, feature] = (clean[feature] - clean[feature].mean()) / std

        formula = "benefit_B ~ " + " + ".join(QC_FEATURES)
        n_projects = clean["project_id"].nunique()
        if n_projects > 1:
            try:
                model = smf.mixedlm(formula, clean, groups=clean["project_id"], re_formula="1")
                result = model.fit(reml=True)
                print(result.summary())
                for name, coef in result.params.items():
                    coeff_rows.append({
                        "model": "mixedlm_benefit_B",
                        "term": name,
                        "coef": coef,
                        "pvalue": result.pvalues.get(name, np.nan),
                    })
            except Exception as exc:
                print(f"MixedLM for benefit_B failed: {exc}")
        else:
            X = sm.add_constant(clean[QC_FEATURES].values)
            y = clean["benefit_B"].values
            result = sm.OLS(y, X).fit()
            terms = ["const"] + QC_FEATURES
            for idx, name in enumerate(terms):
                coeff_rows.append({
                    "model": "ols_benefit_B",
                    "term": name,
                    "coef": result.params[idx],
                    "pvalue": result.pvalues[idx],
                })

    clean_net = df.dropna(subset=QC_FEATURES + ["benefit_B_net"])
    print(f"Rows with complete QC + benefit_B_net data: {len(clean_net)}")
    if len(clean_net) >= 10:
        for feature in QC_FEATURES:
            std = clean_net[feature].std()
            if std and std > 0:
                clean_net.loc[:, feature] = (clean_net[feature] - clean_net[feature].mean()) / std
        formula = "benefit_B_net ~ " + " + ".join(QC_FEATURES)
        try:
            model = smf.mixedlm(formula, clean_net, groups=clean_net["project_id"], re_formula="1")
            result = model.fit(reml=True)
            print(result.summary())
            for name, coef in result.params.items():
                coeff_rows.append({
                    "model": "mixedlm_benefit_B_net",
                    "term": name,
                    "coef": coef,
                    "pvalue": result.pvalues.get(name, np.nan),
                })
        except Exception as exc:
            print(f"MixedLM for benefit_B_net failed: {exc}")

    ordinal_rows = []
    clean_ord = df.dropna(subset=QC_FEATURES + ["t_star_ordinal"])
    if len(clean_ord) >= 10 and clean_ord["t_star_ordinal"].nunique() >= 3:
        try:
            design = clean_ord[QC_FEATURES].copy()
            for feature in QC_FEATURES:
                std = design[feature].std()
                if std and std > 0:
                    design.loc[:, feature] = (design[feature] - design[feature].mean()) / std
            project_dummies = pd.get_dummies(clean_ord["project_id"], prefix="project", drop_first=True)
            design = pd.concat([design, project_dummies], axis=1)
            model = OrderedModel(clean_ord["t_star_ordinal"], design, distr="logit")
            result = model.fit(method="bfgs", disp=False)
            for name, coef in result.params.items():
                ordinal_rows.append({
                    "term": name,
                    "coef": coef,
                    "pvalue": result.pvalues.get(name, np.nan),
                })
            print(result.summary())
        except Exception as exc:
            print(f"Ordered model failed: {exc}")

    write_simple_table(os.path.join(OUT_DIR, "qc_model_results.tsv"), coeff_rows)
    write_simple_table(os.path.join(OUT_DIR, "qc_model_optimal_method.tsv"), ordinal_rows)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plots.")
        return

    plot_df = df.dropna(subset=["benefit_B"] + QC_FEATURES)
    if len(plot_df) >= 5:
        ncols = 2
        nrows = int(np.ceil(len(QC_FEATURES) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4 * nrows))
        axes = np.array(axes).reshape(-1)
        for idx, feature in enumerate(QC_FEATURES):
            ax = axes[idx]
            ax.scatter(plot_df[feature], plot_df["benefit_B"], alpha=0.7, s=28, edgecolors="none")
            ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
            ax.set_xlabel(feature)
            ax.set_ylabel("benefit_B")
            ax.set_title(f"Benefit vs {feature}")
            ax.grid(alpha=0.3)
        for idx in range(len(QC_FEATURES), len(axes)):
            axes[idx].axis("off")
        fig.tight_layout()
        plot_file = os.path.join(PLOTS_DIR, "qc_predictors.png")
        fig.savefig(plot_file, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Written: {plot_file}")

    if coeff_rows:
        coeff_plot_rows = [r for r in coeff_rows if r["model"] == "mixedlm_benefit_B_net"]
        if not coeff_plot_rows:
            coeff_plot_rows = [r for r in coeff_rows if r["model"] == "mixedlm_benefit_B"]
        coeff_plot_rows = [r for r in coeff_plot_rows if r["term"] not in {"Intercept", "Group Var", "const"}]
        if coeff_plot_rows:
            coeff_plot_rows.sort(key=lambda r: float(r["coef"]))
            fig, ax = plt.subplots(figsize=(9, max(4, len(coeff_plot_rows) * 0.45)))
            ax.barh(
                [r["term"] for r in coeff_plot_rows],
                [float(r["coef"]) for r in coeff_plot_rows],
                color="#5C8D89",
            )
            ax.axvline(x=0, color="gray", linestyle="--", alpha=0.6)
            ax.set_xlabel("Coefficient")
            ax.set_title("QC Predictor Coefficients")
            ax.grid(axis="x", alpha=0.3)
            fig.tight_layout()
            plot_file = os.path.join(PLOTS_DIR, "qc_coefficients.png")
            fig.savefig(plot_file, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"Written: {plot_file}")

    print("Training predictive RandomForestClassifier for CLI tool...")
    try:
        from sklearn.ensemble import RandomForestClassifier
        import joblib
        
        # We predict the raw string (U, A, P5...) directly
        rf_df = df.dropna(subset=QC_FEATURES + ["t_star"])
        if len(rf_df) > 0:
            X = rf_df[QC_FEATURES].values
            y = rf_df["t_star"].values
            
            rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight="balanced")
            rf.fit(X, y)
            
            model_path = os.path.join(OUT_DIR, "qc_rf_model.joblib")
            joblib.dump({"model": rf, "features": QC_FEATURES}, model_path)
            print(f"RandomForest model saved to {model_path}")
        else:
            print("Not enough data to train RF model.")
    except Exception as exc:
        print(f"Failed to train RandomForestClassifier: {exc}")

if __name__ == "__main__":
    main()
