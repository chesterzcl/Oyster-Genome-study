import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")
# Load bedtools output
df = pd.read_csv("gene_mean_depth.txt", sep="\t", header=None)

# Set the correct column names
df.columns = ["chr", "start", "end", "gene_id", "dot", "strand", "mean_depth"]



# Optional: log-transform the depth for visualization
df["log10_mean_depth"] = np.log10(df["mean_depth"] + 0.1)
df["gene_length"] = df["end"] - df["start"]
df["length_bin"] = pd.cut(df["gene_length"], bins=[0, 500, 2000, 10000, 50000, 1e8],
                          labels=["<500bp", "0.5–2kb", "2–10kb", "10–50kb", ">50kb"])

median_depth = df["mean_depth"].median()
print(f"Median mean gene depth: {median_depth:.2f}")

# Plot violin plot
plt.figure(figsize=(6, 5))
sns.violinplot(x="length_bin", y="log10_mean_depth", data=df)

# sns.violinplot(y=df["log10_mean_depth"], inner="box", linewidth=1, color="lightblue")
plt.xlabel("Gene Length") 
plt.ylabel("log₁₀(Per Gene Mean Read Depth)")
# plt.title("RNA-seq Support per Gene Model")

#Formatting
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()

# Save figure
plt.savefig("gene_model_support_violin.png", dpi=300)
plt.close()

# Optional: top 10 chromosomes
top_chr = df["chr"].value_counts().nlargest(10).index
df_top = df[df["chr"].isin(top_chr)]
# Map scaffold names to simpler Chr1–Chr10 labels
chr_map = {scaffold: f"Chr{i+1}" for i, scaffold in enumerate(top_chr)}
df_top["chr_label"] = df_top["chr"].map(chr_map)

# Plot
plt.figure(figsize=(6, 5))
sns.violinplot(data=df_top, x="chr_label", y="log10_mean_depth", inner="quartile", scale="width")
plt.xticks(rotation=90)
plt.ylabel("log₁₀(Per Gene Mean Read Depth)")
plt.xlabel("Chromosome")
# plt.title("RNA-seq Support per Gene Model by Chromosome")
plt.tight_layout()
plt.savefig("violin_support_by_chr.png", dpi=300)
plt.close()


# Load annotations
eggnog = pd.read_csv("eggnog_annotation.emapper.annotations", sep="\t", comment="#")

# Optional: rename query column for clarity
eggnog.rename(columns={"query": "gene_id"}, inplace=True)
eggnog["gene_id"] = eggnog["gene_id"].str.replace(r"\.t\d+$", "", regex=True)

# Binary flags for each gene
gene_level = eggnog.groupby("gene_id").agg({
    "eggNOG_OGs": lambda x: any(~x.isin(["-", "None"])),
    "GOs":        lambda x: any(~x.isin(["-", "None", "", pd.NA])),
    "KEGG_ko":    lambda x: any(~x.isin(["-", "None", "", pd.NA])),
    "PFAMs":      lambda x: any(~x.isin(["-", "None", "", pd.NA]))
}).reset_index()

# Count how many genes had each annotation type
total_genes = len(set(df["gene_id"]))
counts = {
    "eggNOG OG": gene_level["eggNOG_OGs"].sum(),
    "GO Terms": gene_level["GOs"].sum(),
    "KEGG Pathways": gene_level["KEGG_ko"].sum(),
    "PFAM Domains": gene_level["PFAMs"].sum(),
    "Unannotated": total_genes - gene_level[["eggNOG_OGs", "GOs", "KEGG_ko", "PFAMs"]].any(axis=1).sum()
}

labels = list(counts.keys())
values = list(counts.values())
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "gray"]

plt.figure(figsize=(6, 5))
bars = plt.bar(labels, values, color=colors)
plt.ylabel("Number of Genes")
plt.xlabel("Annotation Database")
# plt.title("Functional Annotation Coverage per Category")
plt.xticks(rotation=30, ha="right")
plt.grid(axis="y", linestyle="--", alpha=0.5)

# Add counts on top
for bar, count in zip(bars, values):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 300, str(count),
             ha='center', va='bottom', fontsize=9)
    percent = (count / total_genes) * 100
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
             f"{percent:.1f}%", ha='center', va='center', fontsize=9, color="white")    

# Set y-axis limits
plt.ylim(0, max(values) * 1.2)   # starts at 0, adds 20% padding above max
plt.tight_layout()
plt.savefig("panel3_functional_annotation_category.png", dpi=300)
plt.close()



#Segmental duplication
# Load your data
df = pd.read_csv("coords_large_sv.tsv", sep="\t")

# Calculate duplication length
df["dup_length"] = df["ref_end"] - df["ref_start"]

# Clean up category names if needed
df["category"] = df["category"].str.replace("Segmental Duplication \\(Inter-chromosomal\\)", "Inter-chromosomal SD", regex=True)
df["category"] = df["category"].str.replace("Segmental Duplication", "Intra-chromosomal SD")

# Order categories manually if desired
cat_order = ["Tandem Duplication", "Intra-chromosomal Duplication", "Inter-chromosomal Duplication"]

# Plot
plt.figure(figsize=(6, 4))
# Make sure your data is loaded into `df`
df['log10_length'] = np.log10(df['ref_len'])  # Optional for better scaling

plt.figure(figsize=(6, 4))
sns.barplot(
    data=df,
    x='category',
    y='ref_len',  # or 'log10_length' for better scaling
    estimator='median',
    ci='sd',
    palette='Set2'
)
plt.ylabel("Median Duplication Length (bp)")
plt.xlabel("")
plt.title("Median Duplication Length by Category")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig("panel5_duplication_length_barplot.png", dpi=300)
plt.close()


plt.figure(figsize=(6, 4))
sns.countplot(
    data=df,
    x='category',
    palette='Set2'
)
plt.ylabel("Number of Duplications")
plt.xlabel("")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig("panel5_duplication_count_barplot.png", dpi=300)
plt.close()


df['recent'] = df['identity'] > 98

prop_df = df.groupby('category')['recent'].mean().reset_index()

plt.figure(figsize=(6, 4))
sns.barplot(
    data=prop_df,
    x='category',
    y='recent',
    palette='Set2'
)
plt.ylabel("Proportion of Recent Duplications (>98%)")
plt.xlabel("")
plt.xticks(rotation=30, ha="right")
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig("panel5_recent_proportion_barplot.png", dpi=300)
plt.close()


# Start plotting
# Calculate overall Pearson correlation
from scipy.stats import pearsonr

# Compute Pearson correlation
r_val, p_val = pearsonr(df["identity"], df["ref_len"])

# Start plot
plt.figure(figsize=(6, 5))

# Scatter + regression
ax = sns.regplot(
    data=df,
    x='identity',
    y='ref_len',
    scatter_kws={'alpha': 0.5, 's': 30,'color':'purple'},
    line_kws={'color': 'darkblue'}
)

# Axis labels
ax.set_xlabel("Sequence Identity (%)")
ax.set_ylabel("Duplication Length (bp)")

# Add correlation annotation (r and p)
plt.text(
    0.05, 0.95,
    f"r = {r_val:.2f}\np = {p_val:.1e}",
    transform=ax.transAxes,
    verticalalignment='top',
    horizontalalignment='left',
    fontsize=11,
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray')
)

# Clean up visual style
ax.grid(False)

# Save and show
plt.tight_layout()
plt.savefig("panel5_sd_length_identity_clean.png", dpi=300)
plt.show()