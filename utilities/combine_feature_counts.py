import pandas as pd
import os
from functools import reduce

# --- Configuration ---
INPUT_DIR = "/shared/rc/malariarna/ECFC/gene_count"
OUTPUT_FILE = os.path.join(INPUT_DIR, "all_groups_unified_counts.tsv")

# Defined ALL 60 accessions across the three groups
ACCESSIONS_TO_INCLUDE = [
    # Overweight Samples (20 total)
    "SRR24071876", "SRR24071877", 
    "SRR24071895", "SRR24071897", "SRR24071898", "SRR24071899", "SRR24071900",
    "SRR24071901", "SRR24071902", "SRR24071903", "SRR24071904", "SRR24071905",
    "SRR24071906", "SRR24071907", "SRR24071908", "SRR24071909", "SRR24071910",
    "SRR24071911", "SRR24071912", "SRR24071913",
    
    # Lean Samples (20 total)
    "SRR24071916", "SRR24071917", "SRR24071918", "SRR24071919", "SRR24071920",
    "SRR24071921", "SRR24071922", "SRR24071923", "SRR24071924", "SRR24071925",
    "SRR24071926", "SRR24071927", "SRR24071928", "SRR24071929", "SRR24071930",
    "SRR24071931", "SRR24071932", "SRR24071933", "SRR24071934", "SRR24071935",

    # GDM Samples (20 total)
    "SRR24071896", "SRR24071914", "SRR24071915",
    "SRR24071878", "SRR24071879", "SRR24071880", "SRR24071881", "SRR24071882",
    "SRR24071883", "SRR24071884", "SRR24071885", "SRR24071886", "SRR24071887",
    "SRR24071888", "SRR24071889", "SRR24071890", "SRR24071891", "SRR24071892",
    "SRR24071893", "SRR24071894"
]

# Define possible suffixes to check for, in order of preference
COUNT_SUFFIXES = ["_counts.txt", "_dedup_counts.txt"]

# --- Main Script ---

print(f"Starting combination of featureCounts files for ALL THREE groups (Total: {len(ACCESSIONS_TO_INCLUDE)})...")
print(f"Input directory: {INPUT_DIR}")

all_dfs = []
samples_processed = []

for sample_id in ACCESSIONS_TO_INCLUDE:
    file_path = None
    file_type = None

    # Check for preferred files first
    for suffix in COUNT_SUFFIXES:
        temp_path = os.path.join(INPUT_DIR, f"{sample_id}{suffix}")
        if os.path.exists(temp_path):
            file_path = temp_path
            file_type = suffix
            break
    
    if not file_path:
        print(f"WARNING: No suitable count file found for {sample_id} ({COUNT_SUFFIXES}). Skipping.")
        continue

    # Read the featureCounts file. We need the Geneid (column 0) and the Raw Count (column 6).
    try:
        # skiprows=1 is used because featureCounts output has a comment line starting with '#'
        df = pd.read_csv(file_path, sep='\t', skiprows=1, index_col=False)

        # Select the Geneid column (index 0) and the count column (index 6)
        count_column_name = df.columns[6]

        sample_df = df.iloc[:, [0, 6]]
        sample_df = sample_df.rename(columns={
            df.columns[0]: 'Geneid',
            count_column_name: sample_id # Clean column header for Galaxy
        })

        all_dfs.append(sample_df)
        samples_processed.append(sample_id)
        # print(f"Loaded {sample_id} using {file_type} successfully.") # Optional detail logging

    except Exception as e:
        print(f"ERROR processing file {file_path}: {e}. Skipping.")
        continue

if not all_dfs:
    print("FATAL: No count files were successfully loaded. Exiting.")
    exit(1)

# Use functools.reduce and pd.merge to merge all DataFrames based on the 'Geneid' column
print(f"\nMerging {len(samples_processed)} count matrices...")
unified_df = reduce(lambda left, right: pd.merge(left, right, on='Geneid', how='outer'), all_dfs)

# --- Save the Unified Matrix ---
print(f"\nSaving the unified matrix to {OUTPUT_FILE}...")
unified_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

print("\n---------------------------------------------------------")
print(f"SUCCESS! Unified count matrix saved to: {OUTPUT_FILE}")
print(f"Total samples merged: {len(samples_processed)}")
print("---------------------------------------------------------")
