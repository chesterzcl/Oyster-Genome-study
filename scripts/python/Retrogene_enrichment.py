import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from scipy.stats import mannwhitneyu
import seaborn as sns

# === Load BED ===
dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/bed/"
bed_file = "retrogenes.bed"
cols = ['chr', 'start', 'end', 'gene_id']
df = pd.read_csv(dir+bed_file, sep="\t",usecols=[0, 1, 2, 3], header=None, names=cols)


# Clean and compute midpoint
df["start"] = pd.to_numeric(df["start"], errors="coerce")
df["end"] = pd.to_numeric(df["end"], errors="coerce")
df.dropna(subset=["start", "end"], inplace=True)
df["mid"] = ((df["start"] + df["end"]) // 2).astype(int)

# === Set bin size (100 kb) ===
bin_size = 100_000

# === Bin and plot ===
chroms = sorted(df["chr"].unique())
n_chroms = len(chroms)

fig, axes = plt.subplots(nrows=(n_chroms + 1) // 2, ncols=2,
                         figsize=(12, 2.5 * ((n_chroms + 1) // 2)),
                         constrained_layout=True)

axes = axes.flatten()

for i, chrom in enumerate(chroms):
    ax = axes[i]
    chr_df = df[df["chr"] == chrom]
    mids = chr_df["mid"].values

    if len(mids) == 0:
        continue

    max_pos = mids.max()
    bins = np.arange(0, max_pos + bin_size, bin_size)
    counts, _ = np.histogram(mids, bins=bins)

    bin_centers = bins[:-1] + bin_size // 2

    ax.plot(bin_centers, counts, lw=1.5)
    ax.set_title(f"{chrom}")
    ax.set_xlabel("Genomic Position (bp)")
    ax.set_ylabel("Count / 100kb")
    ax.set_xlim([0, max_pos + bin_size])

# Hide unused subplots if any
for j in range(i + 1, len(axes)):
    axes[j].axis("off")

fig.suptitle("Retrogene Density per 100 kb Window", fontsize=14)
plt.savefig("retrogene_density_100kb.png", dpi=300)
plt.show()


all_bins = []

for chrom, grp in df.groupby("chr"):
    mids = grp["mid"].values
    if len(mids) == 0:
        continue

    max_pos = mids.max()
    bins = np.arange(0, max_pos + bin_size, bin_size)
    counts, _ = np.histogram(mids, bins=bins)

    all_bins.append(pd.DataFrame({
        "chr": chrom,
        "bin_start": bins[:-1],
        "bin_end": bins[:-1] + bin_size,
        "retrogene_count": counts
    }))

density_df = pd.concat(all_bins, ignore_index=True)

# --------------------------------------------
# 2. filter for windows with ≥ 30 retrogenes
# --------------------------------------------
high_blocks = density_df.query("retrogene_count >= 20").copy()

# --------------------------------------------
# 3. show them on screen
# --------------------------------------------
print("\nBlocks with ≥ 30 retrogenes per 100 kb:")
print(high_blocks.to_string(index=False))

# --------------------------------------------
# 4. (optional) save to file
# --------------------------------------------
high_blocks.to_csv("high_retrogene_blocks.tsv", sep="\t", index=False)