import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import os

# === SETUP ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/snp_analysis")

# === LOAD FILES ===
het = pd.read_csv("per_sample_het.het", sep="\t")
depth = pd.read_csv("per_sample_depth.idepth", sep="\t")

# === MERGE FILES ===
merged = pd.merge(het, depth, on="INDV", suffixes=("_HET", "_DEPTH"))

# === CALCULATE HETEROZYGOSITY RATE ===
merged["Heterozygosity_Rate"] = (merged["N_SITES_HET"] - merged["O(HOM)"]) / merged["N_SITES_HET"]

# === EXTRACT GROUP INFO (Optional) ===
merged["Group"] = merged["INDV"].str.extract(r"(\d+[LS])")  # e.g., 10L or 11S

# === SELECT AND RENAME COLUMNS ===
summary_df = merged[[
    "INDV", "Group", "N_SITES_HET", "O(HOM)", "E(HOM)", "Heterozygosity_Rate", "F", "MEAN_DEPTH"
]].copy()

summary_df.columns = [
    "Sample", "Group", "Total_SNPs", "Observed_Homozygous", "Expected_Homozygous",
    "Heterozygosity_Rate", "Inbreeding_Coefficient", "Mean_Depth"
]

# === SAVE TABLE ===
summary_df.to_csv("panel1_sample_summary_table.tsv", sep="\t", index=False)




# Optional: Preview
print(summary_df.head())

# Load data
snp_df = pd.read_csv("snp_density_200000.tsv", sep="\t")
pi_df = pd.read_csv("pi_200000.windowed.pi", sep="\t")

# Fix 1-based vs 0-based bin start inconsistency
pi_df["BIN_START"] = pi_df["BIN_START"].astype(int) - 1
pi_df["BIN_END"] = pi_df["BIN_END"].astype(int)

# Merge
merged_df = pd.merge(snp_df, pi_df, on=["CHROM", "BIN_START", "BIN_END"])

# Calculate correlation
r_val, p_val = pearsonr(merged_df["N_VARIANTS"], merged_df["PI"])

# Plot
plt.figure(figsize=(6, 5))
ax = sns.regplot(
    data=merged_df,
    x="N_VARIANTS",
    y="PI",
    scatter_kws={'alpha': 0.5},
    line_kws={'color': 'red'}  # regression line
)

# Labels
plt.xlabel("SNP Count per 200 kb")
plt.ylabel("Nucleotide Diversity (π)")

# Annotate correlation and p-value
# Format p-value more readably
if p_val < 1e-300:
    p_text = "< 1e-300"
else:
    p_text = f"{p_val:.1e}"

plt.text(
    0.05, 0.95,
    f"r = {r_val:.2f}\np = {p_text}",
    transform=ax.transAxes,
    ha="left", va="top",
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray")
)

plt.tight_layout()
plt.savefig("panel2_snp_vs_pi_scatter.png", dpi=300)
plt.close()

# Optional: Correlation summary
correlation = merged_df[["N_VARIANTS", "PI"]].corr().iloc[0, 1]
print(f"Correlation between SNP count and π: {correlation:.3f}")

# Basic stats
summary_stats = {
    "Mean SNPs per bin": merged_df["N_VARIANTS"].mean(),
    "SD SNPs per bin": merged_df["N_VARIANTS"].std(),
    "Min SNPs per bin": merged_df["N_VARIANTS"].min(),
    "Max SNPs per bin": merged_df["N_VARIANTS"].max(),
    "Total SNPs (summed)": merged_df["N_VARIANTS"].sum(),
    "Mean pi": merged_df["PI"].mean(),
    "SD pi": merged_df["PI"].std(),
    "Correlation SNPs vs pi": merged_df[["N_VARIANTS", "PI"]].corr().iloc[0, 1]
}

summary_df = pd.DataFrame.from_dict(summary_stats, orient="index", columns=["Value"])
summary_df.to_csv("panel2_snp_pi_summary.tsv", sep="\t")

print(summary_df)

# Optional: categorize SNP density
bins = merged_df["N_VARIANTS"]
merged_df["SNP_DENSITY_CATEGORY"] = pd.cut(
    bins,
    bins=[-1, bins.quantile(0.33), bins.quantile(0.66), bins.max()],
    labels=["Low", "Medium", "High"]
)

# Count number of bins in each category
category_counts = merged_df["SNP_DENSITY_CATEGORY"].value_counts().sort_index()
print("\nSNP density categories:")
print(category_counts)

# Optional: chromosome-wise mean SNP density
chr_summary = merged_df.groupby("CHROM")["N_VARIANTS"].agg(["mean", "std", "min", "max", "count"])
chr_summary.columns = ["Mean_SNPs", "SD_SNPs", "Min_SNPs", "Max_SNPs", "Num_Bins"]
chr_summary.to_csv("panel2_snp_density_per_chr.tsv", sep="\t")

