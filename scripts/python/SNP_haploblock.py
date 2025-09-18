import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === SETUP ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/snp_analysis")

# === Load .blocks.det file ===
block_file = "oyster_haplo_blocks.blocks.det"  # Change path if needed

# Read only the first 6 columns (ignore haplotype string)
df = pd.read_csv(block_file, sep=r"\s+", header=0, usecols=range(6),
                 names=["CHR", "BP1", "BP2", "KB", "NSNPS", "SNPS"])

# Convert data types if needed
df["KB"] = pd.to_numeric(df["KB"], errors="coerce")
df["NSNPS"] = pd.to_numeric(df["NSNPS"], errors="coerce")

# === Summary statistics ===
total_blocks = len(df)
mean_block_size = df["KB"].mean()
median_block_size = df["KB"].median()
max_block_size = df["KB"].max()
mean_snp_per_block = df["NSNPS"].mean()
median_snp_per_block = df["NSNPS"].median()

print("🧬 Haplotype Block Summary:")
print(f"Total blocks: {total_blocks:,}")
print(f"Mean block size: {mean_block_size:.2f} kb")
print(f"Median block size: {median_block_size:.2f} kb")
print(f"Max block size: {max_block_size:.2f} kb")
print(f"Mean SNPs per block: {mean_snp_per_block:.2f}")
print(f"Median SNPs per block: {median_snp_per_block:.2f}")

# === Plot 1: Block size histogram ===
plt.figure(figsize=(6, 4))
sns.histplot(df["KB"], bins=50, color="steelblue")
plt.xlabel("Block size (kb)")
plt.ylabel("Count")
plt.title("Distribution of haplotype block sizes")
plt.tight_layout()
plt.savefig("block_size_histogram.png", dpi=300)
plt.show()

# === Plot 2: SNP count per block ===
plt.figure(figsize=(6, 4))
sns.histplot(df["NSNPS"], bins=30, color="darkorange")
plt.xlabel("Number of SNPs per block")
plt.ylabel("Count")
plt.title("SNP count per haplotype block")
plt.tight_layout()
plt.savefig("snp_per_block_histogram.png", dpi=300)
plt.show()

# === Optional: Boxplot of block sizes by chromosome ===
plt.figure(figsize=(10, 5))
sns.boxplot(data=df, x="CHR", y="KB", showfliers=False)
plt.ylabel("Block size (kb)")
plt.xticks(rotation=90)
plt.title("Haplotype block sizes by chromosome")
plt.tight_layout()
plt.savefig("block_size_per_chr.png", dpi=300)
plt.show()


# === Per-chromosome summary table ===
summary = df.groupby("CHR").agg(
    Total_Blocks=("KB", "count"),
    Mean_Block_Size_kb=("KB", "mean"),
    Median_Block_Size_kb=("KB", "median"),
    Max_Block_Size_kb=("KB", "max"),
    Mean_SNPs_per_Block=("NSNPS", "mean"),
    Median_SNPs_per_Block=("NSNPS", "median")
).reset_index()

# Round to 3 decimals
summary = summary.round(3)

# Save summary table
summary.to_csv("haplotype_block_summary_per_chr.tsv", sep="\t", index=False)





