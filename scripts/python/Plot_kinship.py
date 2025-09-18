import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/snp_analysis")
# Load the relatedness matrix
genome_df = pd.read_csv("cv20_relatedness.genome", delim_whitespace=True)

# Extract unique sample IDs
samples = sorted(set(genome_df['IID1']).union(set(genome_df['IID2'])))

# Extract short IDs (e.g., 10L) and group info (L or S)
def extract_short_id(full_id):
    match = re.match(r"(\d+[LS])_", full_id)
    return match.group(1) if match else full_id

short_ids = {s: extract_short_id(s) for s in samples}
groups = {sid: ('Large' if 'L' in sid else 'Small') for sid in short_ids.values()}

# Create a symmetric matrix of PI_HAT
matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)
for _, row in genome_df.iterrows():
    matrix.loc[row['IID1'], row['IID2']] = row['PI_HAT']
    matrix.loc[row['IID2'], row['IID1']] = row['PI_HAT']
np.fill_diagonal(matrix.values, 1)

# Rename to short ID
matrix.rename(index=short_ids, columns=short_ids, inplace=True)

# Optional: reorder by group
ordered_ids = sorted(matrix.index, key=lambda x: (groups.get(x, ""), x))
matrix = matrix.loc[ordered_ids, ordered_ids]

# Annotation symbols
def get_asterisk(value):
    if value == 1:
        return "***"
    elif 0.5 < value < 1:
        return "**"
    elif 0.25 <= value <= 0.5:
        return "*"
    else:
        return ""

annot = matrix.applymap(get_asterisk)

# Plotting
plt.figure(figsize=(10, 8))
sns.heatmap(matrix, cmap="Reds", vmin=0, vmax=1, linewidths=0.5,
            linecolor='gray', square=True, xticklabels=True, yticklabels=True,
            annot=annot, fmt="", cbar_kws={"shrink": 0.6, "label": "Pi_Hat"})

plt.xticks(rotation=90)
plt.yticks(rotation=0)

# Legend for asterisk symbols
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='white', edgecolor='gray', label='*   : 0.25–0.5'),
    Patch(facecolor='white', edgecolor='gray', label='**  : 0.5–1'),
    Patch(facecolor='white', edgecolor='gray', label='*** : 1 (identical)')
]
plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.25, 1), frameon=False)

plt.title("Kinship matrix for the first 20 samples")
plt.tight_layout()
plt.savefig("cv20_relatedness_grouped_heatmap.png", dpi=300)
plt.show()