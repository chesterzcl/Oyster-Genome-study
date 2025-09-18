import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os

# Set your working dir
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final")

# Define your giant insertion region manually
FLANK_SIZE = 0
INSERTION_CHROM = "HiC_scaffold_1"
INSERTION_START = 5584859-FLANK_SIZE
INSERTION_END   = 18885845+FLANK_SIZE
# Define flanking window size (on each side)


genes_df = pd.read_csv(
    "genes_function.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "gene_id", "gene_name", "strand"]
)

print(f"Loaded {len(genes_df)} genes from BED.")

def classify_region(row):
    if row['chrom'] != INSERTION_CHROM:
        return 'Background'
    if INSERTION_START <= row['start'] <= INSERTION_END:
        return 'Insertion'
    if (INSERTION_START - FLANK_SIZE) <= row['start'] <= INSERTION_START:
        return 'Flank'
    if INSERTION_END <= row['start'] <= (INSERTION_END + FLANK_SIZE):
        return 'Flank'
    return 'Background'

genes_df['region'] = genes_df.apply(classify_region, axis=1)
print(genes_df['region'].value_counts())



ortho_df = pd.read_csv("all_genes_orthogroup_categories.csv")

merged_df = genes_df.merge(
    ortho_df,
    on="gene_id",
    how="left"
)

merged_df['category'] = merged_df['category'].fillna('Unassigned')
merged_df.to_csv("all_genes_labeled_by_region_and_category.csv", index=False)
counts = merged_df.groupby(['region','category']).size().reset_index(name='gene_count')

###========enrichment test=======================================================
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Separate background and insertion
bg_counts = counts[counts["region"] == "Background"].set_index("category")["gene_count"]
lse_counts = counts[counts["region"] == "Insertion"].set_index("category")["gene_count"]

# Compute total genes in each region
total_bg = bg_counts.sum()
total_lse = lse_counts.sum()

results = []
for category in counts["category"].unique():
    a = lse_counts[category]               # in LSE for this category
    b = bg_counts[category]                # in Background for this category
    c = total_lse - a                      # other LSE genes
    d = total_bg - b                       # other background genes

    table = [[a,b],[c,d]]
    odds, p = fisher_exact(table, alternative='greater')  # enrichment

    results.append([category, a, b, odds, p])

# Make results DataFrame
results_df = pd.DataFrame(results, columns=["category","LSE_count","BG_count","odds_ratio","p_value"])

# Multiple testing correction (FDR)
# results_df["FDR"] = multipletests(results_df["p_value"], method='fdr_bh')[1]

print(results_df.sort_values("p_value"))


#####=========================================================================


total_per_region = counts.groupby('region')['gene_count'].sum().reset_index(name='total')
counts_norm = counts.merge(total_per_region, on='region')
counts_norm['proportion'] = counts_norm['gene_count'] / counts_norm['total']

plt.figure(figsize=(10,6))
sns.barplot(data=counts_norm, x='region', y='proportion', hue='category')
plt.ylabel('Proportion')
plt.xlabel('Region')
plt.legend(title='Ortholog Category')
plt.tight_layout()
plt.savefig("ortholog_category_proportions.png")



insertion_genes = merged_df[
    (merged_df['region'] == 'Insertion') &
    (merged_df['gene_name'].notna()) &
    (merged_df['gene_name'] != "")
]['gene_name'].unique()

print(f"Number of unique gene_names in Insertion: {len(insertion_genes)}")

# Save for gProfiler
pd.Series(insertion_genes).to_csv("insertion_region_gene_names.txt", index=False, header=False)



region_start = INSERTION_START - FLANK_SIZE
region_end = INSERTION_END + FLANK_SIZE

region_genes = merged_df[
    (merged_df['chrom'] == INSERTION_CHROM) &
    (merged_df['start'] >= region_start) &
    (merged_df['start'] <= region_end)
].copy()

print(f"Found {len(region_genes)} genes in visualization window.")
region_genes['midpoint'] = (region_genes['start'] + region_genes['end']) // 2


# ================================
# Load background set
# ================================
all_background_genes = merged_df[
    (merged_df['gene_name'].notna()) &
    (merged_df['gene_name'] != "")
]['gene_name'].unique()

print(f"Number of unique gene_names in background: {len(all_background_genes)}")


# ================================
# Run gProfiler enrichment
# ================================
from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)


results = gp.profile(
    organism='hsapiens',
    query=list(insertion_genes),
    user_threshold=0.05,
    significance_threshold_method="fdr",
    background=list(all_background_genes),
    sources=["GO:BP"],
    no_evidences=True
)

# ================================
# Save enrichment results
# ================================
results.to_csv("insertion_region_GO_enrichment.csv", index=False)
print("*** Saved insertion_region_GO_enrichment.csv ***")


# ================================
# Filter significant results
# ================================
sig_results = results[results["significant"] == True]
sig_results.to_csv("insertion_region_GO_enrichment_significant.csv", index=False)
print(f"Found {len(sig_results)} significant GO terms")


sig_results = sig_results.copy()
sig_results["neg_log10_p_value"] = -np.log10(sig_results["p_value"])

# ================================
# Plot top enriched terms
# ================================


if len(sig_results) > 0:

	# Filter significant terms
	plot_df = sig_results.copy()
	plot_df["foreground_freq"] = plot_df["intersection_size"] / plot_df["query_size"]
	plot_df["background_freq"] = plot_df["term_size"] / plot_df["effective_domain_size"]
	plot_df["enrichment_ratio"] = plot_df["foreground_freq"] / plot_df["background_freq"]
	# Sort GO terms by enrichment_ratio
	plot_df = plot_df.sort_values("enrichment_ratio", ascending=False)
	plot_df["name"] = pd.Categorical(plot_df["name"], categories=plot_df["name"], ordered=True)

	plt.figure(figsize=(14, 10))
	sns.scatterplot(
	    data=plot_df,
	    x="enrichment_ratio",
	    y="name",
	    size="intersection_size",
	    hue="neg_log10_p_value",
	    sizes=(50, 400),
	    palette="viridis",
	    edgecolor="black",
	    alpha=0.8
	)
	plt.xlabel("Enrichment ratio")
	plt.ylabel("GO Biological Process")
	
	# Get the legend
	legend = plt.legend(loc='lower right',borderaxespad=0.5,frameon=False)

	# Manually set titles
	for text in legend.get_texts():
	    if 'intersection_size' in text.get_text():
	        text.set_text(text.get_text().replace('intersection_size', 'Gene Count'))
	    if 'neg_log10_p_value' in text.get_text():
	        text.set_text(text.get_text().replace('neg_log10_p_value', '-log10(p-value)'))

	plt.xlim(-1,20)
	plt.tight_layout()
	plt.savefig("GO_enrichment_bubble_plot.png")
else:
    print("No significant enrichment to plot.")



region_genes = region_genes[region_genes['category'] != 'MG_expansion']
# print(region_genes['category'].unique())

custom_order = ["1:1", "CV_specific", "CV_expansion", "Many:Many", "Unassigned"]

# Map categories to numeric positions based on custom_order
category_map = {cat: i for i, cat in enumerate(custom_order)}
region_genes["y_numeric"] = region_genes["category"].map(category_map)

fig, ax = plt.subplots(figsize=(14,4))

sns.scatterplot(
    data=region_genes,
    x='midpoint',
    y='y_numeric',      # numeric positions
    hue='category',      # color by category
    s=15,                # point size
    ax=ax
)

ax.set_title(f"Orthogroup Category Distribution Along Chr{INSERTION_CHROM.split('_')[-1]}")
ax.set_xlabel(f"Position on Chr{INSERTION_CHROM.split('_')[-1]} (Mb)")
ax.axvspan(INSERTION_START, INSERTION_END, color='gray', alpha=0.2, label='CV-specific Region')

# Custom y-tick labels
ax.set_yticks(range(len(custom_order)))
ax.set_yticklabels(["One-to-One", "CV-specific", "CV Expansion", "Many-to-Many", "Unassigned"])
ax.set_ylabel("Orthogroup Category")

# Remove legend
ax.legend([], [], frameon=False)

plt.tight_layout()
plt.savefig("orthogroup_distribution_scatterplot.png", dpi=600)



# Filter to only genes in the insertion region with assigned orthogroup
layout_df = merged_df[
    (merged_df['chrom'] == INSERTION_CHROM) &
    (merged_df['start'] >= INSERTION_START) &
    (merged_df['end'] <= INSERTION_END) &
    (merged_df['orthogroup'].notna())
].copy()

# Sort by start position for plotting
layout_df = layout_df.sort_values("start")

# Assign a color label (or integer) for each orthogroup
layout_df['orthogroup_label'] = layout_df['orthogroup'].astype("category").cat.codes

# Midpoint for plotting
layout_df['midpoint'] = (layout_df['start'] + layout_df['end']) // 2
layout_df['midpoint_mb'] = layout_df['midpoint'] / 1e6
# Optional: count how many genes per orthogroup in this region
og_counts = layout_df['orthogroup'].value_counts()

# Highlight expanded orthogroups
expanded = og_counts[og_counts > 3].index.tolist()  # adjust threshold as needed
layout_df['highlight'] = layout_df['orthogroup'].apply(lambda x: x if x in expanded else 'Other')
layout_df = layout_df[layout_df['highlight'] != 'Other']



# Sort by chromosome and position
# === 1. Prepare data ===
# Ensure sorted by chromosome and position
layout_df = layout_df.sort_values(['chrom', 'midpoint']).copy()

# Only keep highlighted clusters (orthogroups > 3 genes)
cluster_df_list = []

# === 2. Compute stats per cluster ===
for (chrom, og), group in layout_df.groupby(['chrom', 'highlight']):
    # Skip "Other" clusters if any accidentally remain
    if og == "Other":
        continue
    
    group_sorted = group.sort_values('midpoint')

    # Gene count
    num_genes = len(group_sorted)

    # Cluster start/end & span
    cluster_start = group_sorted['midpoint'].min()
    cluster_end = group_sorted['midpoint'].max()
    cluster_span = cluster_end - cluster_start

    # Within-cluster distances
    dists = group_sorted['midpoint'].diff().dropna()
    avg_within = dists.mean() if not dists.empty else np.nan
    median_within = dists.median() if not dists.empty else np.nan
    min_within = dists.min() if not dists.empty else np.nan
    max_within = dists.max() if not dists.empty else np.nan
    std_within = dists.std() if not dists.empty else np.nan
    cv_within = std_within / avg_within if avg_within and avg_within > 0 else np.nan

    # Strand composition
    num_forward = (group_sorted['strand'] == '+').sum()
    num_reverse = (group_sorted['strand'] == '-').sum()
    strand_consistency = max(num_forward, num_reverse) / num_genes

    # Store cluster summary
    cluster_df_list.append({
        'chrom': chrom,
        'orthogroup': og,
        'num_genes': num_genes,
        'cluster_start': cluster_start,
        'cluster_end': cluster_end,
        'cluster_span': cluster_span,
        'avg_within_distance': avg_within,
        'median_within_distance': median_within,
        'min_within_distance': min_within,
        'max_within_distance': max_within,
        'std_within_distance': std_within,
        'cv_within_distance': cv_within,
        'num_forward_strand': num_forward,
        'num_reverse_strand': num_reverse,
        'strand_consistency': strand_consistency
    })

# Convert to DataFrame
cluster_summary = pd.DataFrame(cluster_df_list)

# === 3. Compute between-cluster distances per chromosome ===
cluster_summary = cluster_summary.sort_values(['chrom', 'cluster_start']).reset_index(drop=True)

cluster_summary['prev_cluster_distance'] = np.nan
cluster_summary['next_cluster_distance'] = np.nan

for chrom, sub in cluster_summary.groupby('chrom'):
    idx = sub.index
    cluster_summary.loc[idx, 'prev_cluster_distance'] = sub['cluster_start'] - sub['cluster_end'].shift(1)
    cluster_summary.loc[idx, 'next_cluster_distance'] = sub['cluster_start'].shift(-1) - sub['cluster_end']

from sklearn.cluster import KMeans

# === 4. K-means classification on CV ===
cv_values = cluster_summary[['cv_within_distance']].values

kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
cluster_summary['kmeans_cluster'] = kmeans.fit_predict(cv_values)

# Determine which cluster is "dispersed" (higher CV) vs "clustered" (lower CV)
centers = kmeans.cluster_centers_.flatten()
low_cluster, high_cluster = np.argsort(centers)

cluster_summary['dispersion_class'] = cluster_summary['kmeans_cluster'].map(
    {low_cluster: 'clustered', high_cluster: 'dispersed'}
)
cluster_summary.to_csv("orthogroup_cluster_lse.csv",index=False)

# Preserve original highlight order as it appears in the plot
unique_orthogroups = list(layout_df['highlight'].unique())


# Create mapping: orthogroup → most frequent gene name, capitalized
rep_gene_names = (
    layout_df[layout_df['gene_name'].notna() & (layout_df['gene_name'] != "")]
    .groupby('highlight')['gene_name']
    .agg(lambda x: x.value_counts().idxmax().upper())  # Capitalize
    .reindex(unique_orthogroups)
)

# Drop NaNs from gene names (but preserve y-axis alignment)
valid_labels = ~rep_gene_names.isna()
plot_orthogroups = np.array(unique_orthogroups)[valid_labels]
plot_gene_names = rep_gene_names.dropna().values

# Plot
fig, ax1 = plt.subplots(figsize=(14, 8))

sns.scatterplot(
    data=layout_df,
    x='midpoint_mb',
    y='highlight',
    hue='highlight',
    palette='tab20',
    legend=False,
    s=60,
    ax=ax1
)

ax1.set_xlabel(f"Position on chr{INSERTION_CHROM.split('_')[-1]} (Mb)",fontsize=14)
ax1.set_ylabel("Orthogroup",fontsize=14)

# Add dashed horizontal lines for labeled orthogroups
for og in plot_orthogroups:
    y_pos = unique_orthogroups.index(og)
    ax1.axhline(y=y_pos, linestyle='--', color='gray', linewidth=0.8, alpha=0.6)

# Set x-axis ticks every 1 Mb
x_min = int(layout_df['midpoint_mb'].min())
x_max = int(layout_df['midpoint_mb'].max()) + 1
x_ticks = np.arange(x_min, x_max + 1, 1)
ax1.set_xticks(x_ticks)
ax1.set_xticklabels([f"{x}" for x in x_ticks],fontsize=10)
ax1.tick_params(axis='y', labelsize=10)  # set y-axis tick label font size to 10

# Add right y-axis with representative gene names (only for valid ones)
ax2 = ax1.twinx()
ax2.set_ylim(ax1.get_ylim())
yticks = [unique_orthogroups.index(og) for og in plot_orthogroups]
ax2.set_yticks(yticks)
ax2.set_yticklabels(plot_gene_names,fontsize=10)
ax2.set_ylabel("Representative Gene Name",fontsize=14)
# Remove top and right spines (borders), keep only x and y axis
for spine in ['top', 'right']:
    ax1.spines[spine].set_visible(False)
    ax2.spines[spine].set_visible(False)
plt.tight_layout()
plt.savefig("insertion_orthogroup_layout_annotated.png",dpi=600)




from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

print("\n*** Analyzing orthogroup enrichment in insertion region ***")

# -------------------------
# Filter only CV species and assigned orthogroup
# -------------------------
all_cv_genes = merged_df[
    (merged_df["species"] == "CV") &
    (merged_df["orthogroup"].notna())
]

insertion_cv_genes = merged_df[
    (merged_df["species"] == "CV") &
    (merged_df["orthogroup"].notna()) &
    (merged_df["region"] == "Insertion")
]

print(f"Total CV genes with orthogroup in genome: {len(all_cv_genes)}")
print(f"Total CV genes with orthogroup in insertion: {len(insertion_cv_genes)}")

# -------------------------
# Count orthogroups
# -------------------------
all_counts = all_cv_genes["orthogroup"].value_counts().rename("all_count")
insertion_counts = insertion_cv_genes["orthogroup"].value_counts().rename("insertion_count")

counts_df = pd.concat([all_counts, insertion_counts], axis=1).fillna(0)
counts_df["all_count"] = counts_df["all_count"].astype(int)
counts_df["insertion_count"] = counts_df["insertion_count"].astype(int)

# -------------------------
# Fisher's exact test
# -------------------------
results = []
total_all = all_counts.sum()
total_insertion = insertion_counts.sum()

for og, row in counts_df.iterrows():
    a = row["insertion_count"]
    b = total_insertion - a
    c = row["all_count"] - a
    d = (total_all - total_insertion) - c
    contingency = [[a, b], [c, d]]
    
    if a > 0:
        try:
            odds_ratio, p_value = fisher_exact(contingency, alternative='greater')
        except:
            odds_ratio, p_value = np.nan, 1.0
        results.append((og, a, row["all_count"], odds_ratio, p_value))

enrichment_df = pd.DataFrame(
    results,
    columns=["orthogroup", "insertion_count", "all_count", "odds_ratio", "p_value"]
)

# FDR correction
enrichment_df["p_adj"] = multipletests(enrichment_df["p_value"], method='fdr_bh')[1]
enrichment_df = enrichment_df.sort_values(by="p_adj")

# -------------------------
# Save table
# -------------------------
enrichment_df.to_csv("orthogroup_enrichment_in_insertion.csv", index=False)
print("*** Saved orthogroup_enrichment_in_insertion.csv ***")

# -------------------------
# Plot top enriched orthogroups
# -------------------------
sig_hits = enrichment_df.query("p_adj < 0.05").sort_values(by="insertion_count", ascending=False).head(20)

if len(sig_hits) > 0:
    plt.figure(figsize=(10,6))
    sns.barplot(
        data=sig_hits,
        y='orthogroup',
        x='insertion_count',
        color="steelblue"
    )
    plt.xlabel("Gene count in insertion")
    plt.ylabel("Orthogroup")
    plt.title("Top Enriched Orthogroups in Insertion Region")
    plt.tight_layout()
    plt.savefig("top_enriched_orthogroups_in_insertion.png")
    plt.close()
    print("*** Saved plot: top_enriched_orthogroups_in_insertion.png ***")
else:
    print("No significant orthogroup enrichment detected to plot.")





