import pandas as pd
import glob
import os

# 1. Load Metadata using the same robust logic
def get_clean_meta(file_path):
    df_raw = pd.read_csv(file_path, header=None)
    header_idx = next(i for i, row in df_raw.iterrows() if 'run' in str(row.values).lower())
    df = pd.read_csv(file_path, skiprows=header_idx)
    df.columns = [str(c).strip().lower() for c in df.columns]
    
    # Identify key columns
    id_col = next(c for c in df.columns if 'run' in c)
    sex_col = next(c for c in df.columns if 'sex' in c)
    
    clean = df[[id_col, sex_col]].copy()
    clean.columns = ['SampleID', 'Sex']
    clean['SampleID'] = clean['SampleID'].astype(str).str.strip()
    return clean

meta = get_clean_meta('METADATA_For_Sample_Sex.csv')
meta_ids = set(meta['SampleID'].unique())

# 2. Get IDs from your folder
folder_files = glob.glob("*_counts.txt")
folder_ids = set([os.path.basename(f).replace('_counts.txt', '') for f in folder_files if "_dedup_" not in f])

# 3. Perform the Audit
missing_in_folder = meta_ids - folder_ids
extra_in_folder = folder_ids - meta_ids

print("--- METADATA AUDIT ---")
print(f"Total Unique IDs in Metadata: {len(meta_ids)}")
print(f"Total Count Files in Folder:  {len(folder_ids)}")

if missing_in_folder:
    print("\n[!] IDs in Metadata but NO file found in folder:")
    for m in sorted(missing_in_folder):
        print(f"  - {m}")

if extra_in_folder:
    print("\n[!] Files in Folder but NOT found in Metadata IDs:")
    for e in sorted(extra_in_folder):
        # Let's check if it's a case sensitivity or hidden character issue
        matches = [m for m in meta_ids if e.lower() in m.lower()]
        print(f"  - {e} (Potential near-matches in meta: {matches})")

# 4. Check for 'NaN' or empty Sex values
empty_sex = meta[meta['Sex'].isna() | (meta['Sex'].astype(str).str.lower() == 'nan')]
if not empty_sex.empty:
    print("\n[!] Samples in Metadata with MISSING Sex information:")
    print(empty_sex)
