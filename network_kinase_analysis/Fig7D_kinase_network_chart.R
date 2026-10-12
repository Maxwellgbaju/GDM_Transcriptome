import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# 1. Configuration
PATH = '/shared/rc/malariarna/ECFC/gene_count/MALE_Children_Results/Overweight_vs_Lean'
# Use the WebGestalt Result file here
INPUT_FILE = os.path.join(PATH, 'enrichment_results_wg_result1772553862.txt')
OUTPUT_BASE = os.path.join(PATH, 'webgestalt_male_overweight_vs_lean_350dpi')

FDR_THRESHOLD = 0.1
FONT_SIZE = 9

plt.rcParams.update({'font.size': FONT_SIZE, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

# 2. Load Data
df = pd.read_csv(INPUT_FILE, sep='\t')
df_sorted = df.sort_values(by='normalizedEnrichmentScore', ascending=True)

# 3. Colors
def assign_colors(row):
    nes, fdr = row['normalizedEnrichmentScore'], row['FDR']
    if nes > 0:
        return '#F08080' if fdr < FDR_THRESHOLD else '#FFC0CB'
    return '#87CEEB' if fdr < FDR_THRESHOLD else '#E0FFFF'

df_sorted['color'] = df_sorted.apply(assign_colors, axis=1)

# 4. Plot
fig, ax = plt.subplots(figsize=(12, 10))
ax.barh(range(len(df_sorted)), df_sorted['normalizedEnrichmentScore'], color=df_sorted['color'], edgecolor='black', linewidth=0.5)

# 5. Labels (Diverging Logic)
for i, (name, nes, fdr) in enumerate(zip(df_sorted['description'], df_sorted['normalizedEnrichmentScore'], df_sorted['FDR'])):
    if nes > 0:
        ax.text(-0.1, i, name, ha='right', va='center')
    else:
        ax.text(0.1, i, name, ha='left', va='center')

# 6. Formatting
ax.axvline(0, color='black', linewidth=1.5)
ax.set_yticks([])
ax.set_xlabel('Normalized Enrichment Score (NES)', fontweight='bold')
ax.set_title('WebGestalt GSEA: Pathway Enrichment (Male Overweight vs Lean)', fontweight='bold', pad=20)
ax.xaxis.grid(True, linestyle='--', alpha=0.3)
for spine in ['top', 'right', 'left']: ax.spines[spine].set_visible(False)

plt.tight_layout()
plt.savefig(f"{OUTPUT_BASE}.png", dpi=350)
plt.savefig(f"{OUTPUT_BASE}.eps", format='eps')
print(f"Figures saved in {PATH}")
