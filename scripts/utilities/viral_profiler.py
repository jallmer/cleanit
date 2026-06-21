import sys
import subprocess
import time

counts_file = sys.argv[1]
bam_file = sys.argv[2]
out_file = sys.argv[3]
top_n = int(sys.argv[4])

print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Executing Strict Viral Mathematical Profiling for Top {top_n} Abundances...")

# 1. Load precise Viral Sequence mapping bounds natively
virus_set = set()
try:
    with open('/pc2/users/o/omiks001/hpc-prf-omiks/ja/genomes/virus_annotation.gtf', 'r') as vf:
        for line in vf:
            if line.strip():
                virus_set.add(line.split('\t')[0])
except Exception as e:
    print(f"ERROR loading viral mapping vectors: {e}")
    sys.exit(1)

# 2. Parse STAR ReadsPerGene.out.tab quantification natively
viral_counts = []
try:
    with open(counts_file, 'r') as f:
        for line in f:
            if line.startswith('N_unmapped') or line.startswith('N_multimapping') or line.startswith('N_noFeature') or line.startswith('N_ambiguous'):
                continue
            parts = line.strip().split()
            if not parts:
                continue
            
            gene_id = parts[0]
            # Column 2 specifies Unstranded Read counts natively
            
            if gene_id in virus_set:
                count = int(parts[1])
                viral_counts.append((gene_id, count))
except Exception as e:
    print(f"ERROR calculating STAR counts matrices: {e}")
    sys.exit(1)

# 3. Sort aggressively by Top Abundance computation
viral_counts.sort(key=lambda x: x[1], reverse=True)
top_viruses = viral_counts[:top_n]

print(f"Identified {len(viral_counts)} active viral matches. Isolating strictly top {top_n} computationally...")

# 4. Process Sequence Covering and File Mapping Flatly
with open(out_file, 'w') as f:
    f.write("Rank\tVirusID\tAbundance_Reads\tMean_Depth\tCoverage_Perc\n")
    for idx, (vid, count) in enumerate(top_viruses):
        if count == 0:
            f.write(f"{idx+1}\t{vid}\t0\t0.0\t0.0\n")
            continue
            
        print(f"  Calculating dynamic SLURM depth map array natively for {vid} ...")
        try:
            res = subprocess.check_output(f"samtools coverage -r '{vid}' '{bam_file}'", shell=True).decode('utf-8')
            lines = res.strip().split('\n')
            if len(lines) > 1:
                data = lines[1].split()
                # Default Samtools coverage columns: #rname(0) startpos(1) endpos(2) numreads(3) covbases(4) coverage(5) meandepth(6)
                coverage_perc = data[5]
                mean_depth = data[6]
            else:
                coverage_perc = "0.0"
                mean_depth = "0.0"
        except Exception as e:
            print(f"  ERROR executing samtools struct depth securely: {e}")
            coverage_perc = "ERROR"
            mean_depth = "ERROR"
            
        try:
            cov_val = float(coverage_perc)
        except ValueError:
            cov_val = 0.0

        if cov_val >= 70.0:
            f.write(f"{idx+1}\t{vid}\t{count}\t{mean_depth}\t{coverage_perc}\n")
        else:
            print(f"  Skipping {vid} due to low coverage ({coverage_perc}%)")

print(f"Viral mapping fully successfully cleanly calculated. Exported safely to: {out_file}")
