#!/usr/bin/env python3
import base64
import os

EXTRACTOR_PATH = "/pc2/users/o/omiks001/scripts/analysis/extract_fastqc_stats.py"
MODEL_PATH = "/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis/qc_rf_model.joblib"
OUT_PATH = "/pc2/users/o/omiks001/scripts/analysis/predict_trimming.py"

def main():
    # 1. Base64 Encode the Model
    with open(MODEL_PATH, "rb") as f:
        model_bytes = f.read()
    model_b64 = base64.b64encode(model_bytes).decode("ascii")

    # 2. Extract parsing logic from extract_fastqc_stats.py
    with open(EXTRACTOR_PATH, "r") as f:
        extractor_code = f.read()

    # Extract lines 41 to 201 (0-indexed 40 to 201)
    lines = extractor_code.split("\n")
    # find def parse_range_width
    start_idx = 0
    for i, line in enumerate(lines):
        if line.startswith("def parse_range_width"):
            start_idx = i
            break
            
    end_idx = 0
    for i in range(start_idx, len(lines)):
        if line.startswith("def weighted_merge"):
            end_idx = i
            break
    
    # Just hardcode the slice based on known lines
    # It starts around def parse_range_width and ends before def weighted_merge
    parse_functions_block = []
    in_block = False
    for line in lines:
        if line.startswith("def parse_range_width"):
            in_block = True
        if line.startswith("def weighted_merge"):
            in_block = False
        if in_block:
            parse_functions_block.append(line)
            
    parsing_code = "\n".join(parse_functions_block)

    # 3. Build the new script
    new_script = f'''#!/usr/bin/env python3
"""
predict_trimming.py

A standalone CLI tool to predict the optimal RNA-Seq trimming strategy 
using a pre-trained RandomForest model from Leave-One-Sample-Out biological concordance data.

This file is a standalone package. It contains the embedded FastQC parser 
and the base64-encoded Random Forest model.
"""

import sys
import os
import pathlib
import warnings
import json
import argparse
import base64
import io
from collections import defaultdict
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import joblib

# ==============================================================================
# EMBEDDED FASTQC PARSER
# ==============================================================================
{parsing_code}

# ==============================================================================
# EMBEDDED MODEL BINARY
# ==============================================================================
MODEL_B64 = "{model_b64}"

def load_embedded_model():
    model_bytes = base64.b64decode(MODEL_B64)
    buffer = io.BytesIO(model_bytes)
    return joblib.load(buffer)

def main():
    help_text = """
Description:
  Predict the mathematically optimal RNA-Seq trimming stringency directly from 
  raw FastQC data using a Random Forest model trained on biological concordance.

Usage Examples:
  1. Standard Pipeline Mode (Outputs pure JSON):
     python3 predict_trimming.py /path/to/fastqc_data.txt
     
  2. Human-Readable Mode (Outputs detailed explanation):
     python3 predict_trimming.py /path/to/fastqc_data.txt --explain
    """
    parser = argparse.ArgumentParser(
        description="Predict optimal RNA-Seq trimming strategy from FastQC data.",
        epilog=help_text,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("fastqc_file", type=str, help="Path to the fastqc_data.txt file.")
    parser.add_argument("--explain", action="store_true", help="Print intricate details and interpretation of the analysis instead of JSON.")
    args = parser.parse_args()

    fastqc_file = pathlib.Path(args.fastqc_file)
    if not fastqc_file.exists():
        if args.explain:
            print(f"Error: File '{{fastqc_file}}' not found.")
        else:
            print(json.dumps({{"error": f"File '{{fastqc_file}}' not found."}}))
        sys.exit(1)

    if args.explain:
        print(f"Parsing FastQC report: {{fastqc_file.name}} ...")
        
    try:
        parsed_data = parse_fastqc_data(fastqc_file)
    except Exception as e:
        if args.explain:
            print(f"Failed to parse FastQC data: {{e}}")
        else:
            print(json.dumps({{"error": f"Failed to parse FastQC data: {{e}}"}}))
        sys.exit(1)

    if args.explain:
        print("Loading embedded predictive QC model...")
        
    model_data = load_embedded_model()
    rf = model_data["model"]
    features = model_data["features"]

    # Convert to DataFrame
    df = pd.DataFrame([parsed_data])
    
    # Handle missing values by coercing to numeric and filling with 0
    for f in features:
        if f not in df.columns:
            df[f] = 0.0
        df[f] = pd.to_numeric(df[f], errors="coerce").fillna(0.0)

    X = df[features].values

    # Predict
    prediction = str(rf.predict(X)[0])
    probs = rf.predict_proba(X)[0]
    
    # Sort probabilities descending
    class_probs = list(zip(rf.classes_, probs))
    class_probs.sort(key=lambda x: x[1], reverse=True)

    if not args.explain:
        # Strict JSON output
        out_dict = {{
            "prediction": prediction,
            "confidence": {{str(c): float(p) for c, p in class_probs}}
        }}
        print(json.dumps(out_dict, indent=2))
    else:
        # Verbose Human-readable output
        print("=====================================================")
        print(f" PREDICTED OPTIMAL TRIMMING STRINGENCY:  {{prediction}} ")
        print("=====================================================")
        print("Confidence Distribution:")
        for c, p in class_probs:
            print(f"  {{c:>4}}:  {{p:>6.2%}}")
            
        print("\\nInterpretation:")
        if prediction == "U":
            print("The model predicts that this sample is of high quality or that trimming ")
            print("will provide little to no biological benefit and risks introducing alignment bias. ")
            print("Recommendation: Proceed with Untrimmed (U) reads.")
        elif prediction == "NA":
            print("The model predicts that this sample's quality profile is highly anomalous or poor.")
            print("It has a high likelihood of catastrophic failure in downstream biological concordance.")
            print("Recommendation: Investigate raw data manually before proceeding.")
        else:
            print(f"The model detected QC artifacts that may negatively impact downstream pathway ")
            print(f"concordance. Trimming using method '{{prediction}}' is recommended to ")
            print("rescue biological signal.")

if __name__ == "__main__":
    main()
'''

    with open(OUT_PATH, "w") as f:
        f.write(new_script)

    print(f"Successfully packaged {OUT_PATH}")

if __name__ == "__main__":
    main()
