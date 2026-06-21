import joblib
import pandas as pd
import numpy as np
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestClassifier

def main():
    # Load the features from the newly built table
    df = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/sample_feature_table.tsv', sep='\t')

    # The model artifact tells us which features were used
    artifact = joblib.load('/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/qc_rf_model.joblib')
    QC_FEATURES = artifact['features']

    # Drop missing values
    df = df.dropna(subset=QC_FEATURES + ['t_star', 'SRR_ID']).copy()

    X = df[QC_FEATURES].values
    y = df['t_star'].values

    # Recreate the exact same RF architecture to get cross-val predictions
    # This ensures we don't overestimate accuracy by testing on the training data
    rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
    y_pred_cv = cross_val_predict(rf, X, y, cv=5)

    # Save predictions
    res = pd.DataFrame({
        'SRR_ID': df['SRR_ID'],
        'True_t_star': y,
        'Predicted_t_star': y_pred_cv
    })
    
    out_path = '/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/concordance/qc_model_predictions.tsv'
    res.to_csv(out_path, sep='\t', index=False)
    print(f"Saved {len(res)} CV predictions to {out_path}")

if __name__ == "__main__":
    main()
