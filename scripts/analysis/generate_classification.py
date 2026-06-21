import pandas as pd
import numpy as np
from pathlib import Path

# Paths
ANALYSIS = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis")
CONCORDANCE = ANALYSIS / "concordance"
IN_FILE = CONCORDANCE / "bio_concordance.tsv"
OUT_FILE = ANALYSIS / "trimming_classification.tsv"

df = pd.read_csv(IN_FILE, sep="\t")

# Pivot to get jaccard_deg and jaccard_pathway for each method per SRR
pivoted_deg = df.pivot_table(index=["project", "SRR_ID"], columns="method", values="jaccard_deg").reset_index()
pivoted_path = df.pivot_table(index=["project", "SRR_ID"], columns="method", values="jaccard_pathway").reset_index()

pivoted = pd.merge(pivoted_deg, pivoted_path, on=["project", "SRR_ID"], suffixes=('_deg', '_path'))

# If 'U_deg' is missing, we can't classify
pivoted = pivoted.dropna(subset=["U_deg", "U_path"])

methods = ["U", "A", "P5", "P10", "P20", "P35"]
THRESHOLD = 0.01

res = []
for _, row in pivoted.iterrows():
    srr_res = {
        "SRR_ID": row["SRR_ID"],
        "project_id": row["project"]
    }
    
    # Determine t_star
    # A method is considered if it exists for both deg and path
    vals_deg = {m: row[f"{m}_deg"] for m in methods if pd.notnull(row.get(f"{m}_deg"))}
    vals_path = {m: row[f"{m}_path"] for m in methods if pd.notnull(row.get(f"{m}_path"))}
    
    if not vals_deg or not vals_path:
        continue
        
    u_deg = row["U_deg"]
    u_path = row["U_path"]
    
    # Determine t_star for Gene-only
    t_star_gene = "U"
    max_val_deg = u_deg
    for m in methods:
        if m == "U" or m not in vals_deg:
            continue
        if vals_deg[m] - u_deg > THRESHOLD and vals_deg[m] > max_val_deg:
            t_star_gene = m
            max_val_deg = vals_deg[m]
            
    # Determine t_star for Pathway-only
    t_star_path = "U"
    max_val_path = u_path
    for m in methods:
        if m == "U" or m not in vals_path:
            continue
        if vals_path[m] - u_path > THRESHOLD and vals_path[m] > max_val_path:
            t_star_path = m
            max_val_path = vals_path[m]
            
    # Determine t_star for Strict BOTH (Option 2)
    t_star_both = "U"
    max_val_both = u_deg
    for m in methods:
        if m == "U" or m not in vals_deg or m not in vals_path:
            continue
        delta_deg = vals_deg[m] - u_deg
        delta_path = vals_path[m] - u_path
        
        if delta_deg > THRESHOLD and delta_path > THRESHOLD:
            if vals_deg[m] > max_val_both:
                t_star_both = m
                max_val_both = vals_deg[m]
    
    # Global t_star (using both for utmost strictness)
    t_star = t_star_both
            
    srr_res["t_star"] = t_star
    srr_res["t_star_gene"] = t_star_gene
    srr_res["t_star_path"] = t_star_path
    
    # Classify each method for gene, path, and both
    for m in methods:
        if m not in vals_deg or m not in vals_path:
            continue
            
        delta_deg = vals_deg[m] - u_deg
        delta_path = vals_path[m] - u_path
        
        # Gene classification
        if m == "U":
            cls_gene = "neutral"
        elif delta_deg > THRESHOLD:
            cls_gene = "helpful"
        elif delta_deg < -THRESHOLD:
            cls_gene = "harmful"
        else:
            cls_gene = "neutral"
            
        # Pathway classification
        if m == "U":
            cls_path = "neutral"
        elif delta_path > THRESHOLD:
            cls_path = "helpful"
        elif delta_path < -THRESHOLD:
            cls_path = "harmful"
        else:
            cls_path = "neutral"
            
        # Both classification (for t_star alignment)
        if m == "U":
            cls_both = "neutral"
        elif delta_deg > THRESHOLD and delta_path > THRESHOLD:
            cls_both = "helpful"
        elif delta_deg < -THRESHOLD and delta_path < -THRESHOLD:
            cls_both = "harmful"
        else:
            cls_both = "neutral"
            
        srr_res[f"{m}_class_gene"] = cls_gene
        srr_res[f"{m}_class_path"] = cls_path
        srr_res[f"{m}_class"] = cls_both # Used for backwards compatibility
        
        srr_res[f"{m}_delta_path"] = delta_deg # keeping backwards compatible name
        srr_res[f"{m}_read_loss"] = np.nan
        srr_res[f"{m}_gene_loss"] = np.nan
        srr_res[f"{m}_assigned_frac_loss"] = np.nan

    res.append(srr_res)

df_out = pd.DataFrame(res)
df_out.to_csv(OUT_FILE, sep="\t", index=False)
print(f"Generated classification for {len(df_out)} samples and saved to {OUT_FILE}")
