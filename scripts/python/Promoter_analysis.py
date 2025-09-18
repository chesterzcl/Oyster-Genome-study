import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Load file ===
dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/"
file_path = "atac_peaks_promoters_annotated.tsv"
df = pd.read_csv(dir+file_path, sep="\t", header=None,
                 names=["Chr", "Peak_Start", "Peak_End", "Support", "Gene_ID", "Strand"])

# === Calculate peak length ===
df["Peak_Length"] = df["Peak_End"]-df["Peak_Start"]

# === Count number of peaks per gene ===
gene_peak_counts = df["Gene_ID"].value_counts().reset_index()
gene_peak_counts.columns = ["Gene_ID", "Num_Peaks"]

# === Summary statistics ===
print("Total unique genes with peaks:", gene_peak_counts.shape[0])
print("Top 10 genes with most peaks:")
print(gene_peak_counts.head(10))

# === Plot distribution of number of peaks per gene ===
plt.figure(figsize=(10, 6))
sns.histplot(gene_peak_counts["Num_Peaks"], bins=30, kde=False)
plt.xlabel("Number of ATAC Peaks per Gene Promoter")
plt.ylabel("Number of Genes")
plt.title("Distribution of ATAC Peaks in Promoter Regions")
plt.tight_layout()
plt.show()

# === Optional: Save summary to file ===
gene_peak_counts.to_csv("gene_promoter_peak_counts.tsv", sep="\t", index=False)