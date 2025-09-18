import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
import re

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster")

ip_file = "20cv_df3_hoh_biSNP_filtered_maf05_ann_fe_sig5_sam8.stats"
df = pd.read_csv(ip_file, sep='\t')

# Compute -log10(p)
df['-logp'] = -np.log10(df['p_value'])

# Index for plotting
df['ind'] = range(len(df))

# Extract numeric chromosome/scaffold index
def extract_chr_num(chrom):
    match = re.search(r'(\d+)$', str(chrom))
    return int(match.group(1)) if match else float('inf')  # Put non-numeric at end

df['chr_num'] = df['#Chromosome'].apply(extract_chr_num)

# Sort by numeric chromosome order and position
df = df.sort_values(by=['chr_num', 'position'])

# Save high significance SNPs by threshold
for i in range(2, 7):
    df[df['-logp'] > i][["#Chromosome", "position", "ref_allele", "alt_allele", "annotation", "-logp"]].to_csv(
        f'high_diff_snp_{i}.txt', sep='\t', index=False
    )

# Regroup properly
df_grouped = df.groupby("#Chromosome")
print(df_grouped.size())

fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
colors = ['lightblue', 'grey']
x_labels = []
x_labels_pos = []

# Group after sorting by chr_num and position
df_grouped = df.groupby("#Chromosome", sort=False)

# Plot each chromosome/scaffold group
for num, (name, group) in enumerate(df_grouped):
    ax.scatter(group['ind'], group['-logp'], color=colors[num % len(colors)], s=3)
    mid_pos = group['ind'].iloc[0] + (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2
    x_labels.append(name.replace("HiC_scaffold_", "chr"))  # Optional renaming
    x_labels_pos.append(mid_pos)

# Bonferroni-corrected threshold (adjust the denominator if necessary)
# plt.axhline(y=-np.log10(0.05 / len(df)), color='red', linestyle='dashed', label='Bonferroni cutoff')
plt.axhline(y=5, color='red', linestyle='dashed')


# X-axis labels
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels, rotation=45, ha='right')

# Axes limits and labels
ax.set_xlim([0, len(df)])
ax.set_ylim([0, 10])
ax.set_xlabel('Chromosome')
ax.set_ylabel('-log(p-value)')
plt.title("Manhattan plot for oyster size")

# Legend
plt.legend(loc='upper right')


# Mark top 20 SNPs (red points and non-overlapping labels)
top_snps = df.nlargest(91, '-logp')

# Plot red dots for top SNPs
ax.scatter(top_snps['ind'], top_snps['-logp'], color='red', s=10, zorder=3)

# Save and show
plt.tight_layout()
plt.savefig("CV_size_Manhattan.png", dpi=800, bbox_inches='tight')
plt.show()