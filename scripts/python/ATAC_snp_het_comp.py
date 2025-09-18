import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# === Set working directory ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/SNP_het_analysis")

# === Load BED files with labels ===
def load_bed(file, label):
    df = pd.read_csv(file, sep="\t", header=None, names=["chr", "start", "end", "het"])
    df["group"] = label
    return df

df_list = [
    load_bed("het_in_ATAC_high.bed", "ATAC High Variance"),
    load_bed("het_in_ATAC_low.bed", "ATAC Low Variance"),
    load_bed("het_in_ATAC_high_promoter.bed", "Promoter + ATAC High"),
    load_bed("het_in_ATAC_low_promoter.bed", "Promoter + ATAC Low"),
    load_bed("het_in_promoters2kb.bed", "Promoters ±2kb"),
    load_bed("het_in_genes.bed", "Gene Body"),
    load_bed("site_het.bed", "Genome-wide")
]

df_all = pd.concat(df_list)

# === Plot ===
plt.figure(figsize=(13, 6))
sns.violinplot(data=df_all, x="group", y="het", inner="box", palette="Set2", cut=0)
# sns.stripplot(data=df_all, x="group", y="het", color="black", alpha=0.1, jitter=0.25)

plt.ylabel("Observed Heterozygosity", fontsize=12)
plt.xlabel("")
plt.title("SNP Heterozygosity Across Genomic Features", fontsize=14)
plt.xticks(rotation=30, ha="right")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("snp_het_violin_comparison.png", dpi=300)
plt.show()