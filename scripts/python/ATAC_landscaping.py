import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import scipy.stats as stats

# === Set working directory ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# Matrix of counts
matrix = pd.read_csv("consensus_counts_matrix.txt", sep="\t")


for i in range(10):
	sample_id=matrix.columns[i+3].split('_')[0]
	matrix.rename(columns={matrix.columns[i+3]: sample_id}, inplace=True)
# print(matrix.columns[3:])

# Annotation
annot = pd.read_csv("ATAC_peaks_annotated_with_labels.txt", sep="\t")

merged = matrix.merge(
    annot[["chr", "start", "end", "location", "label"]],
    on=["chr", "start", "end"],
    how="left"
)

sample_cols = list(merged)[3:-2]
sample_cols = [s for s in sample_cols if s != '5L']

# === Step 1: CPM normalization ===
cpm = merged[sample_cols].div(merged[sample_cols].sum(axis=0), axis=1) * 1e6

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

# === Step 3: Remove top 1% high-signal peaks ===
mean_signal = log2.mean(axis=1)
signal_threshold = mean_signal.quantile(0.99)
log2_filtered = log2[mean_signal < signal_threshold]

plt.figure(figsize=(12, 5))
# sns.boxplot(
#     data=log2_filtered,
#     flierprops=dict(marker='o', markersize=2, markerfacecolor='black', markeredgecolor='black')
# )
sns.violinplot(data=log2_filtered, cut=0, scale='width', inner='quartile')
# plt.title("Boxplot: Log2(CPM + 1)")
plt.ylabel("Log2 CPM")
plt.xlabel("Sample ID")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("step2_log2_violinplot.png")
plt.close()

# 4. Keep corresponding annotation columns
filtered_annotations = merged.loc[mean_signal < signal_threshold, ["chr", "start", "end", "location", "label"]].reset_index(drop=True)
# 5. Combine with CPM
filtered_signal = log2_filtered

# 6. Combine correctly
matrix_filtered = pd.concat([filtered_annotations.reset_index(drop=True),
                             filtered_signal.reset_index(drop=True)], axis=1)
# Reorder columns: genomic coords, then samples, then annotations
reordered_cols = ["chr", "start", "end"] + sample_cols + ["location", "label"]
matrix_filtered = matrix_filtered[reordered_cols]
# print(matrix_filtered)

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

# --------------------------------------------------------
# === Function to plot by peak type===
# --------------------------------------------------------
counts = matrix_filtered['location'].value_counts()

# Optional: reorder if you want to control slice order
counts = counts.reindex(['promoter', 'distal', 'gene_body'])

# Extract clean lists for plotting
counts_values = counts.values.tolist()
counts_labels = [f"{label} ({count})" for label, count in zip(counts.index, counts_values)]

# Define colors
colors = ['#fc8d62','#66c2a5', '#8da0cb']

# Plot
plt.figure(figsize=(6,6))
plt.pie(
    counts_values,
    labels=counts_labels,
    autopct='%1.1f%%',
    startangle=140,
    colors=colors,
    wedgeprops={'edgecolor': 'black'}
)
# plt.title('Distribution of Consensus Peaks by Genomic Location')
plt.tight_layout()
plt.savefig('Consensus_Peak_Location_PieChart.png', dpi=300)
# plt.show()
# --------------------------------------------------------
# === Function to analyze variance in one category ===
# --------------------------------------------------------
xlim=3

def analyze_variance(df, category_name, sample_cols,xlim=5e4):
    print(f"\nAnalyzing: {category_name} peaks (n = {len(df)})")
    
    # Compute variance across samples
    df["variance"] = df[sample_cols].var(axis=1)
    
    # Plot histogram of variance
    plt.figure(figsize=(8,5))
    bin_edges = np.linspace(0, xlim, 1000)
    sns.histplot(df["variance"], bins=bin_edges, color="gray")
    plt.title(f"Variance Distribution for {category_name} Peaks")
    plt.xlabel("Variance across samples")
    plt.ylabel("Number of peaks")
    plt.tight_layout()
    plt.xlim(0,xlim)
    plt.savefig(f"Variance_distribution_{category_name}.png")
    plt.close()
    
    # Print quantiles
    print(df["variance"].describe(percentiles=[0.25,0.5,0.75,0.9,0.95,0.99]))
    
    # Save annotated matrix with variance
    df.to_csv(f"Matrix_with_variance_{category_name}.txt", sep="\t", index=False)
    
    return df

promoter_df=matrix_filtered[matrix_filtered["location"]=="promoter"].copy()
distal_df=matrix_filtered[matrix_filtered["location"]=="distal"].copy()
gene_df=matrix_filtered[matrix_filtered["location"]=="gene_body"].copy()

analyze_variance(promoter_df,"promoter",sample_cols,xlim)
analyze_variance(distal_df,"distal",sample_cols,xlim)
analyze_variance(gene_df,"gene_body",sample_cols,xlim)

promoter_var=promoter_df['variance']
distal_var=distal_df['variance']
genebody_var=gene_df['variance']

from scipy.stats import ks_2samp, mannwhitneyu

# Kolmogorov–Smirnov tests
print("KS test Promoter vs Distal:", ks_2samp(promoter_var, distal_var))
print("KS test Promoter vs Gene body:", ks_2samp(promoter_var, genebody_var))
print("KS test Distal vs Gene body:", ks_2samp(distal_var, genebody_var))

# Mann–Whitney U tests (medians)
print("Mann-Whitney Promoter vs Distal:", mannwhitneyu(promoter_var, distal_var))
print("Mann-Whitney Promoter vs Gene body:", mannwhitneyu(promoter_var, genebody_var))
print("Mann-Whitney Distal vs Gene body:", mannwhitneyu(distal_var, genebody_var))

##Var_dist
plt.figure(figsize=(6, 6))
bin_edges = np.linspace(0, xlim, 500)
sns.histplot(promoter_df["variance"], bins=bin_edges,alpha=0.5, color="#FF6F61", label="Promoter", stat='density')
sns.histplot(distal_df["variance"], bins=bin_edges, alpha=0.5,color="#74C476", label="Distal", stat='density')
sns.histplot(gene_df["variance"], bins=bin_edges, alpha=0.5,color="#6BAED6", label="Gene body", stat='density')
plt.legend()
plt.xlim(0,xlim)
# plt.title("Variance distribution across categories")
plt.xlabel("Per-peak Variance",fontsize=13)
plt.ylabel("Density",fontsize=13)
plt.tight_layout()
plt.savefig(f"Variance_distribution_all.png")


##Mean_dist
mean_df = matrix_filtered.copy()
signal_cols = [col for col in mean_df.columns if col not in ['chr','start','end','location','label']]
# print(mean_df)
# Melt for seaborn
melted = mean_df.melt(id_vars=['location'], value_vars=signal_cols,
                 var_name='Sample', value_name='ATAC_signal')

from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
from itertools import combinations

# First: standardize the 'location' labels (strip spaces, lower-case if needed)
melted['location'] = melted['location'].str.strip()

# See what categories you really have
print("Detected categories in data:", melted['location'].unique())

# Get only existing categories
categories = list(melted['location'].unique())

# Automatically create all valid pairs
pairs = list(combinations(categories, 2))
print("Using pairs:", pairs)

# Optional: Perform manual stats for each pair in console
for (cat1, cat2) in pairs:
    vals1 = melted.loc[melted['location'] == cat1, 'ATAC_signal']
    vals2 = melted.loc[melted['location'] == cat2, 'ATAC_signal']
    u, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
    print(f"{cat1} vs {cat2} Mann-Whitney U p-value: {p:.4e}")

# Plot
plt.figure(figsize=(6,6))
ax = sns.boxplot(
    data=melted,
    x='location',
    y='ATAC_signal',
    palette='Set2',
    showfliers=False,
    width=0.4
)
ax.set_xticklabels(['Distal','Promoter','Gene Body'])

# Add statistical annotation only for valid pairs
if pairs:
    annotator = Annotator(ax, pairs, data=melted, x='location', y='ATAC_signal')
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
    annotator.apply_and_annotate()
else:
    print("No valid pairs found for statistical testing.")

# Labels and style
# plt.title("ATAC Signal by Location Category")
plt.ylabel("Normalized ATAC signal",fontsize=13)
plt.xlabel("Genomic Category",fontsize=13)
# plt.ylim(0,12)
plt.tight_layout()

plt.savefig("ATAC_Boxplot_by_Location_with_pvals.png", dpi=300)
plt.close()



distal_df["mean_log2"] = distal_df[sample_cols].mean(axis=1)
# Optional: only keep top 5% most accessible
cutoff = distal_df["mean_log2"].quantile(0.95)
high_signal_distal = distal_df[distal_df["mean_log2"] >= cutoff].copy()


# Sort by chromosome and start
high_signal_distal = high_signal_distal.sort_values(["chr", "start"])

# Define max distance to cluster
max_distance = 12500  # 12.5 kb typical threshold

clusters = []
current_cluster = []

last_chr = None
last_end = None

for _, row in high_signal_distal.iterrows():
    if last_chr == row["chr"] and (row["start"] - last_end) <= max_distance:
        current_cluster.append(row)
    else:
        if current_cluster:
            clusters.append(pd.DataFrame(current_cluster))
        current_cluster = [row]
    last_chr = row["chr"]
    last_end = row["end"]

if current_cluster:
    clusters.append(pd.DataFrame(current_cluster))

# print(clusters)

print(f"Identified {len(clusters)} putative clusters of high-signal distal peaks")

# === Assign cluster IDs ===
clustered_list = []
for i, cluster_df in enumerate(clusters, 1):   # start IDs at 1
    cluster_df = cluster_df.copy()
    cluster_df["cluster_id"] = f"cluster_{i}"
    clustered_list.append(cluster_df)

# Combine back into single DataFrame
clustered_peaks = pd.concat(clustered_list, ignore_index=True)

# Define sample columns (adjust to your actual column names)
signal_cols = sample_cols

# Sum signals within clusters
cluster_sums = clustered_peaks.groupby("cluster_id")[signal_cols].sum().reset_index()

# Compute total signal for ranking
cluster_sums["total_signal"] = cluster_sums[signal_cols].sum(axis=1)

# Sort by total signal
cluster_sums = cluster_sums.sort_values("total_signal", ascending=False).reset_index(drop=True)
print(cluster_sums.head())


clustered_peaks.round(2).to_csv("DistalPeaks_with_ClusterID.txt", sep="\t", index=False)
# cluster_sums.round(2).to_csv("DistalCluster_SummedSignal_Ranked.txt", sep="\t", index=False)

from kneed import KneeLocator

cluster_sums["rank"] = np.arange(1, len(cluster_sums) + 1)

# --- Plot rank vs total signal ---

x = cluster_sums["rank"].values
y = cluster_sums["total_signal"].values

knee_locator = KneeLocator(x, y, curve="convex", direction="decreasing")
knee_rank = knee_locator.knee
print(knee_rank)
print(f"Detected knee at rank: {knee_rank}")

knee_rank=260

# Label super-enhancers
cluster_sums["super_enhancer"] = cluster_sums["rank"] <= knee_rank


# --- Plot with knee marked ---
plt.figure(figsize=(7,5))
plt.plot(x, y, marker=".", linestyle="-", label="Total Signal")
if knee_rank:
    plt.axvline(knee_rank, color="red", linestyle="--", label=f"Knee at rank {knee_rank}")
plt.xlabel("Cluster Rank",fontsize=13)
plt.ylabel("Total Signal",fontsize=13)
# plt.title("Rank-Ordered Cluster Signal with Knee Point")
# plt.yscale("log")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Cluster_RankPlot_with_Knee.png")
plt.close()


se_clusters = cluster_sums[cluster_sums["super_enhancer"] == True].copy()

# Make cluster_id to unique labels map
cluster_labels = (
    clustered_peaks.groupby("cluster_id")["label"]
    .apply(lambda x: sorted(set(x)))
    .reset_index()
)
cluster_labels["joined_labels"] = cluster_labels["label"].apply(lambda lst: ";".join(lst))

# Restrict to just super-enhancer clusters
cluster_labels = cluster_labels[cluster_labels["cluster_id"].isin(se_clusters["cluster_id"])]

# Merge
se_clusters = se_clusters.merge(
    cluster_labels[["cluster_id", "joined_labels"]],
    on="cluster_id",
    how="left"
)

# Load GFF file
gff = pd.read_csv("braker_with_name_and_description.gff3", sep="\t", header=None, comment="#",
                   names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

# Keep only 'gene' rows
genes = gff[gff["type"] == "gene"].copy()

# Extract gene_id, Name, Description
genes["gene_id"] = genes["attributes"].str.extract(r"ID=([^;]+)")
genes["gene_name"] = genes["attributes"].str.extract(r"Name=([^;]+)")
genes["description"] = genes["attributes"].str.extract(r"Description=([^;]+)")

# Drop duplicates if any
genes = genes.drop_duplicates(subset=["gene_id"])

# Make a dictionary: gene_id -> (gene_name, description)
gene_info_dict = genes.set_index("gene_id")[["gene_name", "description"]].to_dict(orient="index")

def enrich_labels(label_string):
    if pd.isna(label_string):
        return ""
    
    # Split on semicolon to get regions like 'g11176-g11177'
    regions = label_string.split(";")
    genes_seen = set()
    enriched_parts = []

    for region in regions:
        if not region or pd.isna(region):
            continue

        # Split region on hyphen to get individual gene IDs
        gene_ids = region.split("-")
        for gene_id in gene_ids:
            gene_id = gene_id.strip()
            if gene_id in genes_seen:
                continue
            genes_seen.add(gene_id)

            if gene_id in gene_info_dict:
                info = gene_info_dict[gene_id]
                parts = [gene_id]
                if pd.notna(info["gene_name"]):
                    parts.append(info["gene_name"])
                if pd.notna(info["description"]):
                    parts.append(info["description"])
                enriched_parts.append("[{}]".format(";".join(parts)))
            else:
                enriched_parts.append(f"[{gene_id}]")

    return ";".join(enriched_parts)


se_clusters["gene_annotation"] = se_clusters["joined_labels"].apply(enrich_labels)

# Save
se_clusters.round(2).to_csv("super_enhancers_with_gene_labels.txt", sep="\t", index=False)


se_bed_df = clustered_peaks[
    clustered_peaks["cluster_id"].isin(se_clusters['cluster_id'])
][["chr", "start", "end", "label"]]

# Convert columns
se_bed_df["start"] = se_bed_df["start"].astype(int)
se_bed_df["end"] = se_bed_df["end"].astype(int)

# Write to BED
se_bed_df.to_csv("SE_peaks.bed", sep="\t", index=False, header=False)

# Compute mean and variance across your samples
se_only = se_clusters[se_clusters["super_enhancer"]].copy()
se_only["mean_signal"] = se_only[sample_cols].mean(axis=1)
se_only["variance_signal"] = se_only[sample_cols].var(axis=1)

# print(se_clusters)
# Print summary stats
print(se_only[["mean_signal", "variance_signal"]].describe())

from scipy.stats import pearsonr
# Compute correlation
r_val, p_val = pearsonr(se_only["mean_signal"], se_only["variance_signal"])
print(f"Pearson r = {r_val:.3f}, p-value = {p_val:.3e}")

def extract_gene_name(annotation):
    parts = annotation.split('];[')
    ann=""
    for part in parts:
        fields=part.split(';')
        print(fields)
        if len(fields) >= 3:
            if len(ann)!=0:
                ann+=","+fields[1]
            else:
            	ann+=fields[1]
    return ann

se_only["label"] = se_only["gene_annotation"].apply(extract_gene_name)

top_mean = se_only.nlargest(10, "mean_signal")
top_variance = se_only.nlargest(10, "variance_signal")

# Combine these without duplicates
to_label = pd.concat([top_mean, top_variance]).drop_duplicates()
print(to_label[["label", "mean_signal", "variance_signal"]])

# Plot
plt.figure(figsize=(7,5))
# sns.set(style="whitegrid")
ax = sns.regplot(
    data=se_only,
    x="mean_signal",
    y="variance_signal",
    scatter_kws={'alpha': 0.6,'s':10},
    line_kws={'color':'darkred'},
    ci=95
)

offset = -1
for _, row in to_label.iterrows():
    ax.text(
        row["mean_signal"],
        row["variance_signal"]+offset,
        str(row["label"]).upper(),
        fontsize=6,
        weight='bold',
        alpha=0.85,
        ha='center',
        va='top'
    )

# Annotate correlation on the plot
ax.text(
    0.05, 0.95,
    f"Pearson r = {r_val:.2f}\np = {p_val:.3e}",
    transform=ax.transAxes,
    fontsize=10,
    verticalalignment='top',
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.5)
)

ax.set_xlabel("Mean Accessibility", fontsize=13)
ax.set_ylabel("Per-region Accessiblity Variance", fontsize=13)
# ax.set_title("Mean vs Variance in Super-Accessible Regions", fontsize=14)
plt.tight_layout()
plt.savefig("SECluster_Mean_vs_Variance_Annotated.png", dpi=300)
plt.close()



pearson_r, pearson_p = stats.pearsonr(se_only["mean_signal"], se_only["variance_signal"])
print(f"Pearson r = {pearson_r:.3f}, p-value = {pearson_p:.3e}")
spearman_r, spearman_p = stats.spearmanr(se_only["mean_signal"], se_only["variance_signal"])
print(f"Spearman rho = {spearman_r:.3f}, p-value = {spearman_p:.3e}")



# ---------------------------------------------
# Super-enhancer Length Distribution
# ---------------------------------------------
# Get all peaks in SE clusters with location info
se_all_peaks = clustered_peaks[
    clustered_peaks["cluster_id"].isin(se_clusters["cluster_id"])
][["chr", "start", "end", "location", "label", "cluster_id"]]

se_all_peaks["length"] = se_all_peaks["end"] - se_all_peaks["start"]

cluster_spans = (
    se_all_peaks
    .groupby('cluster_id')
    .agg(cluster_chr=('chr', 'first'),
         cluster_start=('start', 'min'),
         cluster_end=('end', 'max'))
    .reset_index()
)

cluster_spans['cluster_length'] = cluster_spans['cluster_end']-cluster_spans['cluster_start']

cluster_spans.to_csv("SE_clusters.bed",sep="\t",index=False,columns=["cluster_chr", "cluster_start", "cluster_end"])

plt.figure(figsize=(8, 5))
sns.histplot(cluster_spans['cluster_length'], bins=50, color='purple')
# plt.title("Super-accessible regulatory hub span length")
plt.xlabel("Cluster Span (bp)")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("SE_Cluster_Length_Histogram.png")
plt.close()

print("\n[INFO] Saved SE_Length_Histogram.png")
print(cluster_spans['cluster_length'].describe(percentiles=[0.25,0.5,0.75,0.9,0.95,0.99]))



########################################################
# === Annotate full consensus peaks with top5% flag ===
########################################################

print("\n[INFO] Annotating all consensus peaks with 'top5pct' flag...")

# Make unique key for all peaks
merged["peak_id"] = merged["chr"].astype(str) + ":" + merged["start"].astype(int).astype(str) + "-" + merged["end"].astype(int).astype(str)
high_signal_distal["peak_id"] = high_signal_distal["chr"].astype(str) + ":" + high_signal_distal["start"].astype(int).astype(str) + "-" + high_signal_distal["end"].astype(int).astype(str)

# Make set of top5% IDs
top5_set = set(high_signal_distal["peak_id"])
# print(top5_set)

# Add column
merged["top5pct_peak_for_SE"] = merged["peak_id"].apply(lambda x: 1 if x in top5_set else 0)
# print(merged[merged["top5pct_peak_for_SE"]==1])

# Save as consensus_annotated.bed
bed_cols = ["chr", "start", "end", "location", "label", "top5pct_peak_for_SE"]
merged[bed_cols].to_csv("consensus_annotated.bed", sep="\t", index=False, header=False)

print("[INFO] Wrote consensus_annotated.bed with top5pct annotation!")



