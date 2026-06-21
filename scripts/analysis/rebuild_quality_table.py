#!/usr/bin/env python3
"""
Parse raw fastqc_data.txt files to extract detailed QC metrics:
  - adapter_rate: mean adapter content fraction across all positions
  - frac_below_q20/q30: fraction of positions with mean quality < 20/30
  - duplication_rate: total duplication percentage
  - n_content: mean N content across positions
  - gc_deviation: std dev of GC content distribution
  - tail_quality_decay: Q_mean_all - Q_mean_last_10%_positions

This fills the gaps that srr_fastqc_metrics.tsv summary data couldn't provide.
"""

import subprocess, re, os, sys
import pandas as pd
import numpy as np

ANALYSIS = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis"

# 1. Build SRR -> fastqc_data.txt path mapping
print("Finding fastqc_data.txt files...")
r1 = subprocess.run(['find', '/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_fastqc/reports/', 
                      '-name', 'fastqc_data.txt', '-path', '*untrimmed*'], 
                     capture_output=True, text=True)
r2 = subprocess.run(['find', '/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_fastqc/archive_pruned_reports/', 
                      '-name', 'fastqc_data.txt', '-path', '*untrimmed*'], 
                     capture_output=True, text=True)

srr_paths = {}  # SRR -> list of fastqc_data.txt paths (R1/R2)
for line in (r1.stdout + r2.stdout).strip().split('\n'):
    if not line: continue
    m = re.search(r'(SRR\d+)', line)
    if m:
        srr = m.group(1)
        if srr not in srr_paths:
            srr_paths[srr] = []
        srr_paths[srr].append(line)

print(f"Found fastqc_data.txt for {len(srr_paths)} unique SRRs")

# 2. Get classification SRRs to process
cls = pd.read_csv(f"{ANALYSIS}/trimming_classification.tsv", sep="\t")
to_process = set(cls['SRR_ID']) & set(srr_paths.keys())
print(f"Classification SRRs with raw data: {len(to_process)}")

# 3. Parse helper
def parse_fastqc_data(filepath):
    """Parse a single fastqc_data.txt and extract metrics."""
    with open(filepath) as f:
        content = f.read()
    
    result = {}
    
    # Per base sequence quality
    m = re.search(r'>>Per base sequence quality\s+\w+\n#Base\tMean.*?\n(.*?)>>END_MODULE', content, re.DOTALL)
    if m:
        lines = m.group(1).strip().split('\n')
        means = []
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    means.append(float(parts[1]))
                except ValueError:
                    pass
        if means:
            result['q_mean'] = np.mean(means)
            result['q_median'] = np.median(means)
            n_pos = len(means)
            result['frac_below_q20'] = sum(1 for q in means if q < 20) / n_pos
            result['frac_below_q30'] = sum(1 for q in means if q < 30) / n_pos
            # Tail quality decay: mean of all - mean of last 10%
            last_10 = means[max(0, int(n_pos * 0.9)):]
            result['tail_quality_decay'] = np.mean(means) - np.mean(last_10)
    
    # Adapter Content
    m = re.search(r'>>Adapter Content\s+\w+\n#Position.*?\n(.*?)>>END_MODULE', content, re.DOTALL)
    if m:
        lines = m.group(1).strip().split('\n')
        adapter_fracs = []
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    # Sum all adapter columns (there may be multiple adapter types)
                    total = sum(float(x) for x in parts[1:])
                    adapter_fracs.append(total / 100.0)  # Convert from % to fraction
                except ValueError:
                    pass
        if adapter_fracs:
            result['adapter_rate'] = np.mean(adapter_fracs)
    
    # Sequence Duplication Levels
    m = re.search(r'#Total Deduplicated Percentage\t([\d.]+)', content)
    if m:
        dedup_pct = float(m.group(1))
        result['duplication_rate'] = (100.0 - dedup_pct) / 100.0  # fraction duplicated
    
    # Per base N content
    m = re.search(r'>>Per base N content\s+\w+\n#Base\tN-Count\n(.*?)>>END_MODULE', content, re.DOTALL)
    if m:
        lines = m.group(1).strip().split('\n')
        n_vals = []
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    n_vals.append(float(parts[1]))
                except ValueError:
                    pass
        if n_vals:
            result['n_content'] = np.mean(n_vals) / 100.0  # Convert from % to fraction
    
    # Per sequence GC content (for gc_deviation)
    m = re.search(r'>>Per sequence GC content\s+\w+\n#GC Content\tCount\n(.*?)>>END_MODULE', content, re.DOTALL)
    if m:
        lines = m.group(1).strip().split('\n')
        gc_vals = []
        counts = []
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    gc_vals.append(float(parts[0]))
                    counts.append(float(parts[1]))
                except ValueError:
                    pass
        if gc_vals and counts:
            total_count = sum(counts)
            if total_count > 0:
                weighted_mean = sum(g * c for g, c in zip(gc_vals, counts)) / total_count
                weighted_var = sum(c * (g - weighted_mean) ** 2 for g, c in zip(gc_vals, counts)) / total_count
                result['gc_deviation'] = np.sqrt(weighted_var) / 100.0  # Normalize to 0-1 scale
                result['gc_content'] = weighted_mean / 100.0
    
    # Basic stats
    m = re.search(r'Total Sequences\t(\d+)', content)
    if m:
        result['total_sequences'] = int(m.group(1))
    
    m = re.search(r'Sequence length\t(\d+)', content)
    if m:
        result['read_length'] = int(m.group(1))
    
    return result

# 4. Process all SRRs
print("Parsing raw FastQC files...")
rows = []
for i, srr in enumerate(sorted(to_process)):
    if i % 100 == 0:
        print(f"  {i}/{len(to_process)}...")
    
    paths = srr_paths[srr]
    # Parse all R1/R2 files and average
    all_metrics = []
    for p in paths:
        try:
            metrics = parse_fastqc_data(p)
            if metrics:
                all_metrics.append(metrics)
        except Exception as e:
            pass
    
    if not all_metrics:
        continue
    
    # Average across R1/R2
    avg = {}
    all_keys = set()
    for m in all_metrics:
        all_keys.update(m.keys())
    
    for k in all_keys:
        vals = [m[k] for m in all_metrics if k in m]
        if vals:
            if k in ('total_sequences',):
                avg[k] = sum(vals)  # Sum for total sequences
            else:
                avg[k] = np.mean(vals)  # Average for everything else
    
    proj = cls[cls['SRR_ID'] == srr]['project_id'].values[0]
    avg['SRR_ID'] = srr
    avg['project_id'] = proj
    rows.append(avg)

df_new = pd.DataFrame(rows)
print(f"Parsed metrics for {len(df_new)} SRRs")

# 5. Rebuild per_srr_quality.tsv with full coverage
# Map column names
col_map = {
    'q_mean': 'Q_mean',
    'q_median': 'Q_median', 
    'total_sequences': 'sequence_depth',
    'read_length': 'read_length_mean',
}

for old, new in col_map.items():
    if old in df_new.columns:
        df_new = df_new.rename(columns={old: new})

# Standardize column names
df_new = df_new.rename(columns={
    'frac_below_q20': 'frac_below_q20',
    'frac_below_q30': 'frac_below_q30',
    'adapter_rate': 'adapter_rate',
    'duplication_rate': 'duplication_rate',
    'n_content': 'n_content',
    'gc_deviation': 'gc_deviation',
    'gc_content': 'gc_content',
    'tail_quality_decay': 'tail_quality_decay',
})

# For the 112 SRRs without raw data, fall back to srr_fastqc_metrics.tsv
fqc = pd.read_csv('/pc2/users/o/omiks001/hpc-prf-omiks/ja/final_fastqc/srr_fastqc_metrics.tsv', sep='\t')
still_missing = set(cls['SRR_ID']) - set(df_new['SRR_ID'])
print(f"\nSRRs without raw data (fallback to summary): {len(still_missing)}")

fallback_rows = []
for _, row in fqc[fqc['srr_id'].isin(still_missing)].iterrows():
    fallback_rows.append({
        'SRR_ID': row['srr_id'],
        'project_id': row['project_id'],
        'sequence_depth': row['total_sequences'],
        'Q_mean': row['mean_quality_all_bases'],
        'Q_median': row.get('median_quality_terminal_bin', np.nan),
        'gc_content': row['gc_percent'] / 100.0 if pd.notnull(row.get('gc_percent')) else np.nan,
        'tail_quality_decay': row['mean_quality_all_bases'] - row['mean_quality_terminal_bin'] if pd.notnull(row.get('mean_quality_terminal_bin')) else np.nan,
        'read_length_mean': float(str(row.get('sequence_length', '0')).split('-')[0]) if pd.notnull(row.get('sequence_length')) else np.nan,
    })
df_fallback = pd.DataFrame(fallback_rows)

# Combine
df_final = pd.concat([df_new, df_fallback], ignore_index=True)

# Ensure all expected columns exist
expected_cols = ['SRR_ID', 'project_id', 'sequence_depth', 'Q_mean', 'Q_median', 
                 'read_length_mean', 'frac_below_q20', 'frac_below_q30',
                 'tail_quality_decay', 'adapter_rate', 'duplication_rate', 
                 'n_content', 'gc_content', 'gc_deviation']
for c in expected_cols:
    if c not in df_final.columns:
        df_final[c] = np.nan

df_final = df_final[expected_cols]
df_final.to_csv(f"{ANALYSIS}/per_srr_quality.tsv", sep="\t", index=False)

# Report
print(f"\nFinal per_srr_quality.tsv: {len(df_final)} SRRs")
print(f"\nFeature coverage:")
for c in expected_cols[2:]:
    n = df_final[c].notna().sum()
    print(f"  {c:25s} {n}/{len(df_final)} ({100*n/len(df_final):.0f}%)")
