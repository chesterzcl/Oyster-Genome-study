import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Load sparse matrix
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")


matrix_file="ALL_100000_matrix.txt"
chrom_sizes_file ="chrom.sizes"
bin_size = 100000

# === Load chromosome sizes ===
chrom_df = pd.read_csv(chrom_sizes_file, sep="\t", header=None, names=["chrom", "length"])
chrom_df = chrom_df[chrom_df["chrom"].str.match(r"HiC_scaffold_[0-9]+")]
chrom_df = chrom_df.sort_values("chrom", key=lambda x: x.str.extract(r"(\d+)").astype(int)[0])

chrom_df["num_bins"] = (chrom_df["length"] // bin_size).astype(int)
chrom_df["start_bin"] = chrom_df["num_bins"].cumsum().shift(fill_value=0)
chrom_df["end_bin"] = chrom_df["start_bin"] + chrom_df["num_bins"]

# === Build global bin index ===
bin_map = {}  # key: (chrom, bin_start) → global_bin_idx
for row in chrom_df.itertuples(index=False):
    for i in range(row.num_bins):
        local_bin = i * bin_size
        global_bin = int(row.start_bin + i)
        bin_map[(row.chrom, local_bin)] = global_bin

n_bins = chrom_df["num_bins"].sum()
contact_matrix = np.zeros((n_bins, n_bins), dtype=np.float32)

# === Load and populate matrix ===
df = pd.read_csv(matrix_file, sep="\t", header=None,
                 names=["chr1", "chr2", "bin1", "bin2", "count"])

missing = 0
for row in df.itertuples(index=False):
    try:
        i = bin_map[(row.chr1, row.bin1)]
        j = bin_map[(row.chr2, row.bin2)]
        contact_matrix[i, j] = row.count
        contact_matrix[j, i] = row.count  # symmetric
    except KeyError:
        missing += 1
print(f"Skipped {missing} unmapped entries")

# === Plot ===

# Transform: square root
data_to_plot = np.log1p(contact_matrix)**3.0

data_to_plot = data_to_plot / np.max(data_to_plot)

# Clip at 99th percentile (excluding zeros)
vmax = np.percentile(data_to_plot[data_to_plot>0],98)

boundaries = chrom_df["start_bin"].tolist()
midpoints = (chrom_df["start_bin"] + chrom_df["num_bins"] // 2).tolist()
labels = [f"chr{i}" for i in range(1, 11)]

# Plot
fig, ax = plt.subplots(figsize=(12,12))
cax = ax.imshow(data_to_plot, cmap="YlOrRd", origin="lower",vmax=vmax)

# === Custom colorbar ticks and labels
norm_ticks = np.linspace(0,1,5)               
real_ticks = norm_ticks*vmax                 

cb = fig.colorbar(cax,ax=ax,fraction=0.05, pad=0.01, ticks=real_ticks,shrink=0.3)
cb.ax.set_yticklabels(["0","0.25","0.5", "0.75", "1.0"])  
cb.ax.tick_params(labelsize=10)
cb.set_label("Normalized contact intensity", fontsize=8, labelpad=5)

# Add gridlines
for b in boundaries:
    ax.axhline(b, color='black',lw=0.5)
    ax.axvline(b, color='black', lw=0.5)

# Label ticks
ax.set_xticks(midpoints)
ax.set_xticklabels(labels, rotation=0,ha="right", fontsize=11)
ax.set_yticks(midpoints)
ax.set_yticklabels(labels, fontsize=11)

plt.tight_layout()
plt.savefig(f"HiC_contact_map_bin{bin_size}.png", dpi=300)
plt.show()