#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from gprofiler import GProfiler

# === Set working directory ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === Load ATAC peak count matrix ===
counts = pd.read_csv("consensus_counts_matrix.txt", sep="\t", header=0)
print(counts.shape)

sample_cols = list(counts)[3:]
sample_cols = [s for s in sample_cols if s != '5L_st_filtered']

# === Step 1: CPM normalization ===
cpm = counts[sample_cols].div(counts[sample_cols].sum(axis=0), axis=1) * 1e6

# --- Plot CPM distribution ---
plt.figure(figsize=(10, 5))
sns.boxplot(data=cpm)
plt.title("Boxplot: CPM")
plt.ylabel("CPM")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("step1_CPM_boxplot.png")
plt.close()

# === Step 2: Log2 transform ===
log2 = np.log2(cpm + 1)

plt.figure(figsize=(10, 5))
sns.boxplot(data=log2)
plt.title("Boxplot: Log2(CPM + 1)")
plt.ylabel("Log2 CPM")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("step2_log2_boxplot.png")
plt.close()

# === Step 3: Remove top 1% high-signal peaks ===
mean_signal = log2.mean(axis=1)
signal_threshold = mean_signal.quantile(0.99)
log2_filtered = log2[mean_signal < signal_threshold]
counts_filtered = counts.loc[mean_signal < signal_threshold].reset_index(drop=True)

plt.figure(figsize=(8, 5))
sns.histplot(mean_signal, bins=100, kde=True)
plt.axvline(signal_threshold, color='red', linestyle='--', label='99th percentile')
plt.title("Distribution of Peak Mean Signal")
plt.xlabel("Mean Log2 CPM")
plt.ylabel("Peak Count")
plt.legend()
plt.tight_layout()
plt.savefig("step3_mean_signal_filtering.png")
plt.close()

# === Step 4: Add coordinates back ===
log2_filtered[["chr", "start", "end"]] = counts_filtered[["chr", "start", "end"]]

# === Step 5: Compute variance ===
log2_filtered["variance"] = log2_filtered[sample_cols].var(axis=1)


# === Step 6: Classify top 20% as high-variance ===
threshold_80 = log2_filtered["variance"].quantile(0.80)
log2_filtered["variance_class"] = np.where(
    log2_filtered["variance"] >= threshold_80,
    "high",
    "low"
)

print(f"80th percentile variance threshold: {threshold_80:.4f}")

# === Plot variance histogram with threshold ===
plt.figure(figsize=(8, 5))
plt.hist(log2_filtered["variance"], bins=100, color="gray")
plt.axvline(threshold_80, color="red", linestyle="--", label='80th percentile cutoff')
plt.xlabel("Variance across samples")
plt.ylabel("Number of peaks")
plt.title("Variance distribution of ATAC peaks")
plt.legend()
plt.tight_layout()
plt.savefig("Per-peak_variance_distribution_with_threshold.png")
plt.close()

# === Output BED files ===
high_bed = log2_filtered[log2_filtered["variance_class"] == "high"][["chr", "start", "end"]]
low_bed = log2_filtered[log2_filtered["variance_class"] == "low"][["chr", "start", "end"]]

high_bed.to_csv("ATAC_high_variance_peaks.bed", sep="\t", header=False, index=False)
low_bed.to_csv("ATAC_low_variance_peaks.bed", sep="\t", header=False, index=False)

# === Plot sorted variance bar chart ===
sorted_variance = log2_filtered.sort_values("variance").reset_index(drop=True)
plt.figure(figsize=(12, 5))
plt.bar(
    x=np.arange(len(sorted_variance)),
    height=sorted_variance["variance"],
    color=sorted_variance["variance_class"].map({"high": "red", "low": "blue"}),
    width=1.0
)
plt.title("Per-peak Variance Distribution with Top 20% Highlighted")
plt.xlabel("Peaks (sorted by variance)")
plt.ylabel("Variance")
plt.legend(handles=[
    plt.Line2D([0], [0], color='red', lw=4, label='High Variance (Top 20%)'),
    plt.Line2D([0], [0], color='blue', lw=4, label='Low Variance')
])
plt.tight_layout()
plt.savefig("Variance_distribution_quantile_classified.png")
plt.close()

# === Annotate and combine ===
high_var = log2_filtered[log2_filtered["variance_class"] == "high"].copy()
low_var = log2_filtered[log2_filtered["variance_class"] == "low"].copy()
high_var["type"] = "high"
low_var["type"] = "low"
combined_var = pd.concat([high_var, low_var])
combined_var_sorted = combined_var.sort_values("type", ascending=False)

# === Load promoter regions ===
promoters = pd.read_csv("gene_promoters_1k.bed", sep="\t", header=None)
promoters.columns = ["chr", "start", "end", "gene_id", "annotation", "strand"]

def get_gene_id(row, df):
    hits = df[
        (df["chr"] == row["chr"]) &
        (df["end"] >= row["start"]) &
        (df["start"] <= row["end"])
    ]
    return hits["gene_id"].values[0] if not hits.empty else np.nan

high_var["gene_id"] = high_var.apply(lambda r: get_gene_id(r, promoters), axis=1)
low_var["gene_id"] = low_var.apply(lambda r: get_gene_id(r, promoters), axis=1)

high_genes = high_var["gene_id"].dropna().unique()
low_genes = low_var["gene_id"].dropna().unique()

print(f"High-variance gene count: {len(high_genes)}")
print(f"Low-variance gene count: {len(low_genes)}")

# === Map to gene names using GFF ===
gff = pd.read_csv("braker_with_name_and_description.gff3", sep="\t", header=None, comment="#",
                  names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
genes = gff[gff["type"] == "gene"]
genes["gene_id"] = genes["attributes"].str.extract(r"ID=([^;]+)")
genes["gene_name"] = genes["attributes"].str.extract(r"Name=([^;]+)")

gene_map = genes.dropna(subset=["gene_id", "gene_name"]).drop_duplicates(subset="gene_id")
high_named = gene_map[gene_map["gene_id"].isin(high_genes)]["gene_name"].dropna().unique().tolist()
low_named = gene_map[gene_map["gene_id"].isin(low_genes)]["gene_name"].dropna().unique().tolist()

pd.DataFrame({"gene_name": high_named}).to_csv("Genes_high_var.csv", index=False)
pd.DataFrame({"gene_name": low_named}).to_csv("Genes_low_var.csv", index=False)

# === GO enrichment ===
gp = GProfiler(return_dataframe=True)
background_genes = list(set(high_named + low_named))
print(f"High: {len(high_named)}, Low: {len(low_named)}, Background: {len(background_genes)}")

go_high = gp.profile(organism="hsapiens", query=high_named,
                      no_evidences=False, user_threshold=1.0, background=background_genes)
go_low = gp.profile(organism="hsapiens", query=low_named,
                     no_evidences=False, user_threshold=1.0, background=background_genes)

go_high.to_csv("GO_enrichment_high_var.csv", index=False)
go_low.to_csv("GO_enrichment_low_var.csv", index=False)

# === Plot top GO terms ===
top_n = 10

for label, df, color, title in [
    ('high', go_high, 'Reds_r', 'High-Variance Genes'),
    ('low', go_low, 'Blues_r', 'Low-Variance Genes')
]:
    top_terms = df.sort_values("p_value").head(top_n)
    plt.figure(figsize=(8, 6))
    sns.barplot(
        y=top_terms["name"],
        x=-np.log10(top_terms["p_value"]),
        palette=color
    )
    plt.xlabel("-log10(p-value)")
    plt.ylabel("GO Term")
    plt.title(f"Top GO Terms ({title})")
    plt.tight_layout()
    plt.savefig(f"GO_enrichment_{label}_variance.png")
    plt.close()