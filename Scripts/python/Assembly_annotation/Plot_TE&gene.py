import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
import matplotlib.patches as patches


# Parameters
window_size = 100000

# Load inputs
fai = pd.read_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/primary_dedup_chr.fa.fai",sep="\t",header=None,usecols=[0,1],names=["chrom","length"])
genes = pd.read_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/gene_analysis/genes.bed",sep="\t",header=None,names=["chrom","start","end","id"])
tes = pd.read_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/repeat_analysis/TEs.bed",sep="\t",header=None,names=["chrom", "start", "end", "id"])
out_dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/gene_analysis/"

def bin_features(df, fai, window_size, label):
    records = []
    for chrom, length in fai.values:
        sub = df[df["chrom"] == chrom]
        bins = np.arange(0, length + window_size, window_size)
        for start in bins[:-1]:
            end = start + window_size
            overlapping = sub[(sub["end"] > start) & (sub["start"] < end)].copy()
            if not overlapping.empty:
                overlap = overlapping.apply(
                    lambda row: max(0, min(end, row["end"]) - max(start, row["start"])), axis=1
                ).sum()
            else:
                overlap = 0
            records.append([chrom, start, end, overlap / window_size])  # normalized fraction
    return pd.DataFrame(records, columns=["chrom", "start", "end", f"{label}_frac"])

# Bin both TE and gene features
gene_bins = bin_features(genes,fai,window_size, "gene")
te_bins = bin_features(tes,fai,window_size, "te")

# Merge bins
merged = pd.merge(gene_bins,te_bins,on=["chrom", "start", "end"], how="outer").fillna(0)

print(merged)


# Plot
sns.set(style="whitegrid")
plt.figure(figsize=(7, 6))
sns.scatterplot(data=merged, x="te_frac", y="gene_frac", alpha=0.5)
plt.xlabel("TE coverage fraction per window")
plt.ylabel("Gene coverage fraction per window")
plt.title("TE vs Gene Density Relationship")
plt.tight_layout()
plt.savefig(f"{out_dir}te_vs_gene_density.png", dpi=300)
# plt.show()


# Extract TE and gene coverage arrays
te_cov = merged["te_frac"]
gene_cov = merged["gene_frac"]

# Calculate Pearson correlation
r, p_value = pearsonr(te_cov, gene_cov)

print(f"Pearson correlation (r): {r:.4f}")
print(f"p-value: {p_value:.4e}")


# Input: merged DataFrame with 'gene_frac' and 'te_frac'
X = merged[["gene_frac", "te_frac"]].values

# Fit KMeans (try 3 or 4 clusters)
kmeans = KMeans(n_clusters=3,random_state=42)
merged["cluster"] = kmeans.fit_predict(X)

# Optional: name clusters based on centroids
centroids = kmeans.cluster_centers_
for i, (gene_c, te_c) in enumerate(centroids):
    print(f"Cluster {i}: gene={gene_c:.2f}, TE={te_c:.2f}")

domain_labels = {
    0: "Mixed",
    1: "Gene-rich",
    2: "TE-rich"
}
merged["domain"] = merged["cluster"].map(domain_labels)

domain_summary = merged.groupby("domain")["end"].count().reset_index(name="windows")
domain_summary["Mb_covered"] = domain_summary["windows"] * (window_size/1e6)
domain_summary["percent"] = 100 * domain_summary["Mb_covered"] / domain_summary["Mb_covered"].sum()
print(domain_summary)

merged[["chrom","start","end","domain"]].to_csv(f"{out_dir}genome_domains_3class.bed",sep="\t",index=False,header=False)


sns.set(style="whitegrid")
plt.figure(figsize=(7, 6))
sns.scatterplot(data=merged,x="te_frac",y="gene_frac",hue="domain",palette="Set2",alpha=0.6)

# Annotate centroids for each domain
centroids = merged.groupby("domain")[["te_frac","gene_frac"]].mean()
for label, row in centroids.iterrows():
    plt.text(row["te_frac"]+0.01,row["gene_frac"]+0.01, label,
             fontsize=11,weight="bold",bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="gray"))

# Axes and title
plt.xlabel("TE Coverage Fraction (per window)")
plt.ylabel("Gene Coverage Fraction (per window)")
plt.title("Genome Functional Domains by Gene vs TE Density")
plt.legend(title="Domain", loc="best")
plt.tight_layout()
plt.savefig(f"{out_dir}te_gene_domains_labeled.png",dpi=300)





domains=merged[["chrom","start","end","domain"]]

fai['chrom_num'] = fai['chrom'].str.extract(r'_(\d+)$').astype(int)
fai = fai.sort_values('length', ascending=False).reset_index(drop=True)

# === Color mapping for domains ===
# pick a nice 3-color palette
palette = sns.color_palette("Set2", 3)  
# or: palette = plt.get_cmap("tab10").colors[:3]

domain_colors = {
    'Gene-rich': palette[0],    # soft green
    'TE-rich':   palette[1],    # warm orange
    'Mixed':     palette[2]     # gentle purple
}

# === Plot setup ===
fig, ax = plt.subplots(figsize=(12, len(fai)*0.35))
chrom_height = 0.5
# Remove top/left/right spines
for s in ['top','left','right']:
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_linewidth(1)
ax.get_yaxis().set_visible(False)
ax.grid(False)
plt.subplots_adjust(left=0.15, right=0.85, bottom=0.15)

# === Draw chromosomes and domains ===
for idx, row in fai.iterrows():
    chrom = row['chrom']
    length_mb = row['length'] / 1e6
    y = idx

    # background bar
    ax.add_patch(patches.Rectangle((0, y), length_mb, chrom_height,
                                   facecolor='lightgray', edgecolor='black', linewidth=0.5))
    # overlay domain segments
    sub = domains[domains['chrom']==chrom]
    for _, feat in sub.iterrows():
        start_mb = feat['start'] / 1e6
        end_mb   = feat['end']   / 1e6
        ax.add_patch(patches.Rectangle((start_mb, y), end_mb-start_mb, chrom_height,
                                       facecolor=domain_colors.get(feat['domain'],'gray'), edgecolor='none'))
    # chromosome label
    label = f"chr{row['chrom_num']}"
    ax.text(-1.5, y+chrom_height/2, label, va='center', ha='right', fontsize=8)


# force the y-limits to cover every row
ax.set_ylim(-0.5, len(fai) - 0.5)

# now label and save
ax.set_xlim(-2, fai['length'].max() / 1e6 + 2)
ax.set_xlabel('Position (Mb)')
handles = [patches.Patch(color=c, label=k) for k, c in domain_colors.items()]
ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.95,1), frameon=False)

plt.savefig(f"{out_dir}TE_gene_dist_chromosome.png", dpi=300)
plt.show()









# High gene, low TE
gene_rich = merged[(merged["gene_frac"]>0.6)&(merged["te_frac"]<0.1)]

# High TE, low gene
te_rich = merged[(merged["te_frac"]>0.4)&(merged["gene_frac"]<0.1)]

# Dual extremes (both high)
conflict = merged[(merged["gene_frac"]>0.6)&(merged["te_frac"]>0.4)]

# Deserts (both low)
deserts = merged[(merged["gene_frac"]<0.05)&(merged["te_frac"]<0.05)]

gene_rich[["chrom", "start", "end"]].to_csv(f"{out_dir}gene_rich_islands.bed", sep="\t", index=False, header=False)
te_rich[["chrom", "start", "end"]].to_csv(f"{out_dir}te_rich_islands.bed", sep="\t", index=False, header=False)
conflict[["chrom", "start", "end"]].to_csv(f"{out_dir}high_gene_te_regions.bed", sep="\t", index=False, header=False)
deserts[["chrom", "start", "end"]].to_csv(f"{out_dir}genome_deserts.bed", sep="\t", index=False, header=False)

print(f"Gene-rich windows: {len(gene_rich)}")
print(f"TE-rich windows: {len(te_rich)}")
print(f"High-gene & high-TE windows: {len(conflict)}")
print(f"Desert windows: {len(deserts)}")

