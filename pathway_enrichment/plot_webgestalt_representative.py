import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 1. Configuration
# Updated input filename as requested
INPUT_FILE = 'enrichment_results_wg_result1771516143.txt'
OUTPUT_FILE_PNG = 'webgestalt_female_replica_350dpi.png'
OUTPUT_FILE_EPS = 'webgestalt_female_replica_350dpi.eps'
FDR_THRESHOLD = 0.1
FONT_SIZE = 9

# Set global font properties to 9pt Arial
plt.rcParams.update({
    'font.size': FONT_SIZE,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
})

# 2. Load and Prepare Data
try:
    df = pd.read_csv(INPUT_FILE, sep='\t')
except FileNotFoundError:
    print(f"Error: {INPUT_FILE} not found. Please ensure the file is in the script directory.")
    exit()

# Sort by NES so the chart flows from most negative to most positive
df_sorted = df.sort_values(by='normalizedEnrichmentScore', ascending=True)

# 3. Define WebGestalt Color Logic
def assign_colors(row):
    nes = row['normalizedEnrichmentScore']
    fdr = row['FDR']
    if nes > 0:
        # WebGestalt Salmon/Coral for Up-regulated
        return '#F08080' if fdr < FDR_THRESHOLD else '#FFC0CB'
    else:
        # WebGestalt Sky Blue for Down-regulated
        return '#87CEEB' if fdr < FDR_THRESHOLD else '#E0FFFF'

df_sorted['color'] = df_sorted.apply(assign_colors, axis=1)

# 4. Create the Figure
fig, ax = plt.subplots(figsize=(12, 10))

# Create the bars
bars = ax.barh(range(len(df_sorted)), 
               df_sorted['normalizedEnrichmentScore'], 
               color=df_sorted['color'], 
               edgecolor='black', 
               linewidth=0.5,
               height=0.7)

# 5. Position Labels like WebGestalt (Diverging)
for i, (name, nes, fdr) in enumerate(zip(df_sorted['description'], 
                                         df_sorted['normalizedEnrichmentScore'], 
                                         df_sorted['FDR'])):
    # The Pathway Name
    if nes > 0:
        # Text on the left of the 0-line for positive bars
        ax.text(-0.1, i, name, ha='right', va='center', fontsize=FONT_SIZE)
        # Stats on the right of the bar
        ax.text(nes + 0.05, i, f'FDR: {fdr:.3f}', ha='left', va='center', fontsize=FONT_SIZE-1)
    else:
        # Text on the right of the 0-line for negative bars
        ax.text(0.1, i, name, ha='left', va='center', fontsize=FONT_SIZE)
        # Stats on the left of the bar
        ax.text(nes - 0.05, i, f'FDR: {fdr:.3f}', ha='right', va='center', fontsize=FONT_SIZE-1)

# 6. Aesthetic Formatting
ax.axvline(0, color='black', linewidth=1.5) # The central spine
ax.set_yticks([]) # Hide the default Y-axis ticks
ax.set_xlabel('Normalized Enrichment Score (NES)', fontweight='bold', fontsize=FONT_SIZE)

# Updated Title for Female Cohort
ax.set_title('WebGestalt GSEA: Pathway Enrichment (Female GDM vs Lean)', 
             fontweight='bold', pad=20, fontsize=FONT_SIZE+2)

# Set X-axis limits to make room for labels
limit = max(abs(df_sorted['normalizedEnrichmentScore'])) + 2.5
ax.set_xlim(-limit, limit)

# Vertical grid lines like WebGestalt
ax.xaxis.grid(True, linestyle='--', alpha=0.3)
ax.set_axisbelow(True)

# Remove the outer box (spines) for a cleaner look
for spine in ['top', 'right', 'left']:
    ax.spines[spine].set_visible(False)

# 7. Legend
legend_elements = [
    mpatches.Patch(color='#F08080', label=f'Positive NES (FDR < {FDR_THRESHOLD})'),
    mpatches.Patch(color='#87CEEB', label=f'Negative NES (FDR < {FDR_THRESHOLD})'),
    mpatches.Patch(color='#FFC0CB', label='Not Significant (FDR > 0.1)')
]
ax.legend(handles=legend_elements, loc='upper left', frameon=True, fontsize=FONT_SIZE)

plt.tight_layout()

# 8. Save High Resolution in both PNG and EPS formats
plt.savefig(OUTPUT_FILE_PNG, dpi=350)
plt.savefig(OUTPUT_FILE_EPS, format='eps')

print(f"Success! Female GDM replica saved as {OUTPUT_FILE_PNG} and {OUTPUT_FILE_EPS}")
