import pandas as pd
import os

# 1. Load the existing master file
master_file = 'Merged_Counts_Master.csv'
if os.path.exists(master_file):
    master_df = pd.read_csv(master_file)
    print(f"Loaded existing master file with {len(master_df.columns)-1} samples.")
else:
    print("Error: Merged_Counts_Master.csv not found in the current directory.")
    # Exit or handle as needed
    exit()

# 2. Define the range of samples to add
start_id = 24071906
end_id = 24071935

# 3. Iterate and merge
for i in range(start_id, end_id + 1):
    sample_id = f"SRR{i}"
    file_name = f"{sample_id}_counts.txt"
    
    if os.path.exists(file_name):
        print(f"Processing {file_name}...")
        
        # Read the featureCounts file (skip the first line of metadata)
        temp_df = pd.read_csv(file_name, sep='\t', skiprows=1)
        
        # Extract Geneid and the last column (which contains the actual counts)
        # We rename the count column to the Sample ID for clarity
        counts_col = temp_df.columns[-1]
        temp_df = temp_df[['Geneid', counts_col]]
        temp_df.columns = ['Geneid', sample_id]
        
        # Merge with the master dataframe using an outer join
        master_df = pd.merge(master_df, temp_df, on='Geneid', how='outer')
    else:
        print(f"File {file_name} not found. Skipping...")

# 4. Final Cleanup: Fill any missing genes across samples with 0
master_df = master_df.fillna(0)

# 5. Save the updated table
output_name = 'Merged_Counts_Master_Final.csv'
master_df.to_csv(output_name, index=False)

print(f"\nMerge complete!")
print(f"Total samples in final table: {len(master_df.columns) - 1}")
print(f"Saved to: {output_name}")
