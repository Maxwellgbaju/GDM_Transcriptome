import pandas as pd
import glob
import os

def load_metadata(file_path):
    print(f"Reading {file_path}...")
    df_raw = pd.read_csv(file_path, header=None)
    
    header_idx = None
    for i, row in df_raw.iterrows():
        row_values = row.astype(str).str.lower().tolist()
        if any('run' in str(v) for v in row_values) and any('sex' in str(v) for v in row_values):
            header_idx = i
            break
            
    if header_idx is None:
        raise KeyError("Could not find header row with 'Run' and 'sex' columns.")
    
    df = pd.read_csv(file_path, skiprows=header_idx)
    df.columns = [str(c).strip().lower() for c in df.columns]
    
    col_map = {}
    for col in df.columns:
        if 'run' in col: col_map['id'] = col
        if 'group' in col: col_map['cond'] = col
        if 'sex' in col or 'gender' in col: col_map['sex'] = col
            
    clean_df = df[[col_map['id'], col_map['cond'], col_map['sex']]].copy()
    clean_df.columns = ['SampleID', 'Condition', 'Sex']
    for col in clean_df.columns:
        clean_df[col] = clean_df[col].astype(str).str.strip()
    clean_df = clean_df[clean_df['SampleID'] != 'nan']
    return clean_df

# --- Execution ---

try:
    meta = load_metadata('METADATA_For_Sample_Sex.csv')
    
    # FIND ALL COUNT FILES
    # We look for *_counts.txt but EXCLUDE *_dedup_counts.txt
    all_files = glob.glob("*_counts.txt")
    files = [f for f in all_files if "_dedup_" not in f]
    
    if not files:
        print("Error: No files ending in _counts.txt found (checked for files without '_dedup_')")
    else:
        print(f"Merging {len(files)} standard count files...")
        all_counts = []
        for f in files:
            # Extract SRR ID (e.g., SRR24071876_counts.txt -> SRR24071876)
            srr_id = os.path.basename(f).replace('_counts.txt', '')
            
            # featureCounts format: GeneID is col 0, Count is last col
            df_counts = pd.read_csv(f, sep='\t', skiprows=1)
            sample_counts = df_counts.iloc[:, [0, -1]]
            sample_counts.columns = ['GeneID', srr_id]
            sample_counts.set_index('GeneID', inplace=True)
            all_counts.append(sample_counts)

        master_counts = pd.concat(all_counts, axis=1)

        # Split by Sex
        m_ids = meta[meta['Sex'].str.lower().isin(['male', 'm'])]['SampleID'].tolist()
        f_ids = meta[meta['Sex'].str.lower().isin(['female', 'f'])]['SampleID'].tolist()

        m_cols = [s for s in m_ids if s in master_counts.columns]
        f_cols = [s for s in f_ids if s in master_counts.columns]

        # Save results
        master_counts[m_cols].to_csv("counts_MALES.csv")
        master_counts[f_cols].to_csv("counts_FEMALES.csv")
        meta.to_csv("final_metadata_cleaned.csv", index=False)

        print("-" * 30)
        print(f"DONE!")
        print(f"Total files processed: {len(files)}")
        print(f"Created counts_MALES.csv ({len(m_cols)} samples matched)")
        print(f"Created counts_FEMALES.csv ({len(f_cols)} samples matched)")
        
        # Check for missing matches
        if len(m_cols) + len(f_cols) < len(files):
            print(f"Warning: {len(files) - (len(m_cols) + len(f_cols))} files did not match IDs in metadata.")

except Exception as e:
    print(f"An error occurred: {e}")
