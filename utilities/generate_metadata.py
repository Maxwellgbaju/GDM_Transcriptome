import pandas as pd
import os

# --- Configuration ---
# NOTE: This script assumes the output directory is the same as your count files
INPUT_DIR = "/shared/rc/malariarna/ECFC/gene_count"
OUTPUT_FILE = os.path.join(INPUT_DIR, "all_groups_metadata.tsv")

# Define the samples and their corresponding groups (Condition)
sample_groups = {}

# Overweight Samples (20 total)
overweight_accessions = [
    "SRR24071876", "SRR24071877", 
    "SRR24071895", "SRR24071897", "SRR24071898", "SRR24071899", "SRR24071900",
    "SRR24071901", "SRR24071902", "SRR24071903", "SRR24071904", "SRR24071905",
    "SRR24071906", "SRR24071907", "SRR24071908", "SRR24071909", "SRR24071910",
    "SRR24071911", "SRR24071912", "SRR24071913",
]
for acc in overweight_accessions:
    sample_groups[acc] = "Overweight"

# Lean Samples (20 total)
lean_accessions = [
    "SRR24071916", "SRR24071917", "SRR24071918", "SRR24071919", "SRR24071920",
    "SRR24071921", "SRR24071922", "SRR24071923", "SRR24071924", "SRR24071925",
    "SRR24071926", "SRR24071927", "SRR24071928", "SRR24071929", "SRR24071930",
    "SRR24071931", "SRR24071932", "SRR24071933", "SRR24071934", "SRR24071935"
]
for acc in lean_accessions:
    sample_groups[acc] = "Lean"

# GDM Samples (20 total)
gdm_accessions = [
    "SRR24071896", "SRR24071914", "SRR24071915",
    "SRR24071878", "SRR24071879", "SRR24071880", "SRR24071881", "SRR24071882",
    "SRR24071883", "SRR24071884", "SRR24071885", "SRR24071886", "SRR24071887",
    "SRR24071888", "SRR24071889", "SRR24071890", "SRR24071891", "SRR24071892",
    "SRR24071893", "SRR24071894"
]
for acc in gdm_accessions:
    sample_groups[acc] = "GDM"

# Create the DataFrame
# The columns must be: SampleID, Condition
df_meta = pd.DataFrame(list(sample_groups.items()), columns=['SampleID', 'Condition'])

# IMPORTANT: Sorting by SampleID (ascending) ensures the order matches the
# column names in the unified count matrix, which is necessary for DGE tools.
df_meta = df_meta.sort_values(by='SampleID').reset_index(drop=True)

# --- Save the Metadata File ---
print(f"Generating metadata file for 60 samples...")
print(f"Saving to: {OUTPUT_FILE}")

# Use index=False and sep='\t' for Galaxy compatibility
df_meta.to_csv(OUTPUT_FILE, sep='\t', index=False)

print("\n---------------------------------------------------------")
print(f"SUCCESS! Metadata file saved to: {OUTPUT_FILE}")
print("Ready for upload alongside the unified count matrix.")
print("---------------------------------------------------------")
