import os
import sqlite3
import pandas as pd

DB_PATH = "/pc2/users/o/omiks001/srr_queue.db"
ELIGIBILITY_CSV = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/BioProject_Eligibility_Summary.csv"
OUTPUT_CSV = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/Supplementary_Table_BioProjects.csv"
SAMPLE_SHEETS_DIR = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/deseq2_metadata/sample_sheets"

def main():
    # Load Eligibility Info (This has 48 curated projects)
    eligibility_df = pd.read_csv(ELIGIBILITY_CSV)
    
    # Get pipeline stats from DB
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    
    db_stats = {}
    c.execute("""
        SELECT project_id, 
               count(srr_id) as total_pipeline_srrs,
               sum(CASE WHEN status='done' THEN 1 ELSE 0 END) as compute_successful_srrs,
               sum(CASE WHEN fastq_size_bytes > 30000000000 THEN 1 ELSE 0 END) as large_fastqs
        FROM srr_queue 
        GROUP BY project_id
    """)
    for row in c.fetchall():
        db_stats[row["project_id"]] = {
            "Total_Pipeline_SRRs": row["total_pipeline_srrs"],
            "Compute_Successful_SRRs": row["compute_successful_srrs"],
            "Had_Large_FastQs_>30GB": "Yes" if row["large_fastqs"] > 0 else "No"
        }
        
    # Get organism/platform info from a master list or just sample sheets if needed
    # Actually, we can get Organism/Instrument from the old table or just leave it out if we don't have it here.
    # The previous script might have pulled it. Let's load the old Supplementary table to keep that metadata.
    old_supp_df = pd.read_csv(OUTPUT_CSV)
    metadata_map = {}
    for _, row in old_supp_df.iterrows():
        metadata_map[row["BioProject"]] = {
            "Organism": row["Organism"],
            "Instrument_Model": row["Instrument_Model"],
            "Library_Layout": row["Library_Layout"]
        }

    out_rows = []
    
    for _, row in eligibility_df.iterrows():
        project = row["BioProject"]
        
        # Pipeline metrics
        p_stats = db_stats.get(project, {"Total_Pipeline_SRRs": 0, "Compute_Successful_SRRs": 0, "Had_Large_FastQs_>30GB": "No"})
        
        # Curated metadata counts (from the sample sheets)
        sheet_path = os.path.join(SAMPLE_SHEETS_DIR, f"{project}.csv")
        curated_srrs = 0
        if os.path.exists(sheet_path):
            sheet_df = pd.read_csv(sheet_path)
            curated_srrs = len(sheet_df)
            
        m_info = metadata_map.get(project, {"Organism": "Unknown", "Instrument_Model": "Unknown", "Library_Layout": "Unknown"})
            
        out_rows.append({
            "BioProject": project,
            "Organism": m_info["Organism"],
            "Instrument_Model": m_info["Instrument_Model"],
            "Library_Layout": m_info["Library_Layout"],
            "Total_Pipeline_SRRs": p_stats["Total_Pipeline_SRRs"],
            "Compute_Successful_SRRs": p_stats["Compute_Successful_SRRs"],
            "Had_Large_FastQs_>30GB": p_stats["Had_Large_FastQs_>30GB"],
            "Curated_Metadata_SRRs": curated_srrs,
            "Total_Classes": row["Total_Classes"],
            "Classes_with_Replicates": row["Classes_with_Replicates"],
            "Included_in_PyDESeq2": row["Included_in_PyDESeq2"],
            "Eligible_PyDESeq2_SRRs": row["Eligible_SRRs"] if row["Included_in_PyDESeq2"] == "Yes" else 0,
            "Reason_for_Exclusion": row["Reason_for_Exclusion"]
        })
        
    final_df = pd.DataFrame(out_rows)
    final_df.to_csv(OUTPUT_CSV, index=False)
    
    print(f"Generated {OUTPUT_CSV}")
    print(f"Total Projects: {len(final_df)}")
    print(f"Total Pipeline SRRs: {final_df['Total_Pipeline_SRRs'].sum()}")
    print(f"Total Compute Successful: {final_df['Compute_Successful_SRRs'].sum()}")
    print(f"Total Curated Metadata SRRs: {final_df['Curated_Metadata_SRRs'].sum()}")
    print(f"Total PyDESeq2 Eligible SRRs: {final_df['Eligible_PyDESeq2_SRRs'].sum()}")

if __name__ == "__main__":
    main()
