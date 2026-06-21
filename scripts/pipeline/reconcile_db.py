#!/usr/bin/env python3
"""
reconcile_db.py — Reconcile on-disk outputs with the SQLite tracking queue.

Scans the analysis directory for SRRs and inspects precisely how many of the 6
expected featureCounts matrices exist.
- 6/6 modes found: updates DB status to 'done'.
- 1-5/6 modes found: updates DB status to 'pending' to resume pipeline loops.

Usage:
    python3 scripts/analysis/reconcile_db.py --dry-run
    python3 scripts/analysis/reconcile_db.py --apply
"""

import os
import sys
import glob
import re
import sqlite3
import argparse
from collections import defaultdict

# HPC paths
BASE_DIR = "/scratch/hpc-prf-omiks/ja"
DB_PATHS = [
    os.path.expanduser("~/srr_queue.db"),
    os.path.expanduser("~/scripts/srr_queue.db")
]

def find_active_db():
    for p in DB_PATHS:
        if os.path.exists(p):
            return p
    return None

def analyze_disk_completion(base_dir):
    """
    Returns a dict mapping SRR_ID -> dict(project_id, found_modes_count)
    """
    srr_info = defaultdict(lambda: {"project": "unknown", "modes_found": 0})
    srr_pat = re.compile(r"(SRR\d+)")

    counts_dir = os.path.join(base_dir, "flattened_counts")
    if not os.path.isdir(counts_dir):
        return srr_info

    # 6 possible file patterns per SRR
    # untrmd_{SRR_ID}_fC.txt.gz
    # {SRR_ID}_trimmomatic_adapter_fC.txt.gz
    # {SRR_ID}_trimmomatic_P5_fC.txt.gz
    # {SRR_ID}_trimmomatic_P10_fC.txt.gz
    # {SRR_ID}_trimmomatic_P20_fC.txt.gz
    # {SRR_ID}_trimmomatic_P35_fC.txt.gz
    
    for root, dirs, files in os.walk(counts_dir):
        project = os.path.basename(root)
        if not project.startswith("PRJ"):
            continue
            
        for f in files:
            # We only care about the count matrices
            if not f.endswith("_fC.txt.gz"):
                continue
                
            # Extract SRR id
            m = srr_pat.search(f)
            if m:
                srr = m.group(1)
                srr_info[srr]["project"] = project
                srr_info[srr]["modes_found"] += 1

    return srr_info

def main():
    parser = argparse.ArgumentParser(description="Reconcile DB with on-disk completion states")
    parser.add_argument("--dry-run", action="store_true", help="Print changes without modifying DB")
    parser.add_argument("--apply", action="store_true", help="Apply state to the DB")
    args = parser.parse_args()

    if not args.dry_run and not args.apply:
        print("Please specify either --dry-run or --apply")
        sys.exit(1)

    if not os.path.isdir(BASE_DIR):
        print(f"Error: Base dir not found: {BASE_DIR}")
        sys.exit(1)

    db_path = find_active_db()
    if not db_path:
        print("Error: Could not locate srr_queue.db in expected locations.")
        sys.exit(1)

    print(f"Using database: {db_path}")
    print(f"Scanning disk: {BASE_DIR}/flattened_counts")
    
    disk_srrs = analyze_disk_completion(BASE_DIR)
    print(f"Found {len(disk_srrs)} unique SRRs with Count matrices on disk.")

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Get current DB state
    cur.execute("SELECT srr_id, status FROM srr_queue")
    db_state = {row[0]: row[1] for row in cur.fetchall()}
    print(f"Found {len(db_state)} unique SRRs in DB queue.")

    # Classify disks
    to_mark_done = []
    to_mark_pending = []

    for srr, info in disk_srrs.items():
        modes_count = info["modes_found"]
        current_status = db_state.get(srr, "MISSING")
        
        # If it has all 6, it should be 'done'
        if modes_count >= 6:
            if current_status != "done":
                to_mark_done.append((srr, info["project"], current_status, modes_count))
        # If 1-5, it stalled and should be 'pending' so it can resume
        elif modes_count > 0:
            if current_status not in ["success", "completed", "done"]:
                # Also we don't want to override if it's currently successfully 'submitted' or 'todo'
                # But to guarantee completion, we push to 'pending' if it was missing or failed.
                to_mark_pending.append((srr, info["project"], current_status, modes_count))

    print("\n--- Reconciliation Report ---")
    print(f"SRRs on disk fully complete (6/6 modes) but NOT marked 'done' in DB: {len(to_mark_done)}")
    if to_mark_done:
        print(f"  Examples:")
        for s, p, c_stat, cnt in to_mark_done[:5]:
            print(f"    {s} ({cnt}/6 modes, currently: {c_stat})")

    print(f"\nSRRs partially complete (1-5/6 modes) needing resume ('pending'):  {len(to_mark_pending)}")
    if to_mark_pending:
        print(f"  Examples:")
        for s, p, c_stat, cnt in to_mark_pending[:5]:
            print(f"    {s} ({cnt}/6 modes, currently: {c_stat})")

    if not to_mark_done and not to_mark_pending:
        print("\nAll disk SRRs are perfectly aligned with DB state. Nothing to do!")
        conn.close()
        sys.exit(0)

    if args.apply:
        print("\nApplying changes to database...")
        added_count = 0
        
        # Apply DONE
        for srr, project, _, _ in to_mark_done:
            cur.execute(
                "INSERT OR REPLACE INTO srr_queue (srr_id, project_id, status) VALUES (?, ?, 'done')",
                (srr, project)
            )
            added_count += 1
            
        # Apply PENDING
        for srr, project, _, _ in to_mark_pending:
            # We use pending or todo depending on pipeline preference. Usually 'todo' is the start
            cur.execute(
                "INSERT OR REPLACE INTO srr_queue (srr_id, project_id, status) VALUES (?, ?, 'todo')",
                (srr, project)
            )
            added_count += 1
            
        conn.commit()
        print(f"Successfully adjusted {added_count} records in DB.")
    else:
        print("\n[DRY RUN] Run with --apply to commit these changes to the DB.")

    conn.close()

if __name__ == "__main__":
    main()
