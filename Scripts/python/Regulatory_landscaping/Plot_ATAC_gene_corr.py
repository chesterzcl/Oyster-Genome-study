import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import os

# Set working directory if needed
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === Load Data ===
gene_df = pd.read_csv("gene_density_200k.tsv", sep="\t")
atac_df = pd.read_csv("atac_peak_density_200k.tsv", sep="\t")

# === Merge Data ===
merged_df = pd.merge(
    gene_df,
    atac_df[["chrom", "bin_start", "bin_end", "atac_pct"]],
    on=["chrom", "bin_start", "bin_end"]
)

# === Correlation ===
r, p = pearsonr(merged_df["gene_count"], merged_df["atac_pct"])

# === Plot Setup ===
sns.set_theme(style="white")
plt.figure(figsize=(6, 5))
highlight_color = "#E76F51"

# === Scatter and regression ===
sns.regplot(
    data=merged_df,
    x="gene_count",
    y="atac_pct",
    scatter_kws={"s": 20, "alpha": 1, "color": highlight_color},
    line_kws={"color": "#264653", "linewidth": 2},
)

# === Annotate ===
plt.xlabel("Gene Count per 200 kb", fontsize=13)
plt.ylabel("ATAC-seq Peak Coverage (%)", fontsize=13)
# plt.title("Gene Density vs. Chromatin Accessibility", fontsize=14, pad=12)

plt.text(
    0.05,
    0.95,
    f"Pearson r = {r:.2f}\np-value = {p:.2e}",
    transform=plt.gca().transAxes,
    fontsize=12,
    ha="left",
    va="top",
    bbox=dict(facecolor="white", edgecolor="gray", boxstyle="round,pad=0.3")
)

sns.despine()
plt.grid(False)
plt.tight_layout()
plt.savefig("ATAC_gene_corr_200k.png")
plt.show()