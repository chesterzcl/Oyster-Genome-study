import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# =============================
# SETTINGS
# =============================
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

threshold_FE = 1.2
threshold_log10P = 25

files = {
    'Gene-body': 'knownResults_gene.txt',
    'Promoter': 'knownResults_promoter.txt',
    'Distal': 'knownResults_distal.txt'
}

# =============================
# Function to read and process HOMER output
# =============================
def process_homer_file(filename, category):
    print(f"Processing {filename}...")
    df = pd.read_csv(filename, sep='\t')
    df.columns = df.columns.str.strip()

    # Clean motif name
    df['Motif Name'] = df['Motif Name'].str.split('/').str[0]
    df['TF_Family'] = df['Motif Name'].str.extract(r'\((.*?)\)')

    # Parse % columns robustly
    df['% of Target Sequences with Motif'] = df['% of Target Sequences with Motif'].str.rstrip('%')
    df['% of Background Sequences with Motif'] = df['% of Background Sequences with Motif'].str.rstrip('%')
    df['% of Target Sequences with Motif'] = pd.to_numeric(df['% of Target Sequences with Motif'], errors='coerce')
    df['% of Background Sequences with Motif'] = pd.to_numeric(df['% of Background Sequences with Motif'], errors='coerce')

    # Compute Fold Enrichment
    df['FoldEnrichment'] = df['% of Target Sequences with Motif'] / df['% of Background Sequences with Motif']

    # Correct -log10(p) conversion from HOMER's ln(p)
    df['Log P-value'] = pd.to_numeric(df['Log P-value'], errors='coerce')
    df['-log10p'] = -df['Log P-value'] / np.log(10)

    # Drop bad rows
    df = df.dropna(subset=['FoldEnrichment', '-log10p'])

    # Return all data (no per-file filtering)
    df_filtered = df[['Motif Name', 'TF_Family', 'FoldEnrichment', '-log10p']]
    df_filtered['Category'] = category
    return df_filtered

# =============================
# Process all files
# =============================
all_data = []
for cat, file in files.items():
    df = process_homer_file(file, cat)
    all_data.append(df)

all_data = pd.concat(all_data, ignore_index=True)
print("\n✅ Total motif entries:", len(all_data))

# =============================
# Build TF Family mapping
# =============================
motif_family_map = all_data.drop_duplicates('Motif Name').set_index('Motif Name')['TF_Family']

# =============================
# Pivot both metrics
# =============================
FE_df = all_data.pivot_table(
    index='Motif Name',
    columns='Category',
    values='FoldEnrichment',
    fill_value=0
)

log10p_df = all_data.pivot_table(
    index='Motif Name',
    columns='Category',
    values='-log10p',
    fill_value=0
)

# =============================
# Clean inf and NaN
# =============================
FE_df.replace([np.inf, -np.inf], np.nan, inplace=True)
log10p_df.replace([np.inf, -np.inf], np.nan, inplace=True)
FE_df = FE_df.dropna()
log10p_df = log10p_df.dropna()

# =============================
# Filter: keep motifs passing BOTH thresholds in ANY region
# =============================
pass_FE = (FE_df > threshold_FE)
pass_P = (log10p_df > threshold_log10P)
pass_both = pass_FE & pass_P
keep_rows = pass_both.any(axis=1)

print("\nMotifs retained after threshold filtering:", keep_rows.sum())

# Apply filter
FE_df = FE_df[keep_rows]

# =============================
# Add back TF Family for sorting
# =============================
FE_df['TF_Family'] = motif_family_map
FE_df = FE_df.sort_values(['TF_Family'] + list(files.keys()))

# Drop TF_Family for plotting
FE_df = FE_df.reset_index().set_index('Motif Name')
FE_df = FE_df.drop(columns=['TF_Family'])

# =============================
# Plotting
# =============================
plt.figure(figsize=(8, max(6, 0.2*len(FE_df))))
sns.set(style='whitegrid')

ax = sns.heatmap(
    FE_df,
    cmap='vlag',
    center=1,
    annot=False,
    linewidths=0.1,
    cbar_kws={'label': 'Fold Enrichment', 'shrink': 0.3}
)

ax.set_yticks(np.arange(len(FE_df))+0.5)
ax.set_yticklabels(FE_df.index, rotation=0, fontsize=6)

plt.xlabel('Peak Category', fontsize=13)
plt.ylabel('Motif', fontsize=13)
plt.tight_layout()

plt.savefig('Motif_Enrichment_Heatmap.png', dpi=300)
plt.show()

print("\nHeatmap generated successfully!")


