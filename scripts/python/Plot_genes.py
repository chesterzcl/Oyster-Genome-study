import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# Inputs
gff_file = "braker.gff3"
fai_file = "primary_dedup_chr_masked_hp_sealed.fa.fai"


# Parameters
window_size = 100000

# Load scaffold sizes
fai = pd.read_csv(fai_file,sep="\t",header=None,usecols=[0,1],names=["chrom", "length"])
fai = fai[["chrom", "length"]]

# Load GFF3 and filter gene features
gff=pd.read_csv(gff_file, sep="\t", comment="#", header=None,
                  names=["chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"])
genes=gff[gff["feature"] == "gene"]
# Filter to only genes that are protein-coding



mrna = gff[gff["feature"] == "mRNA"].copy()
mrna["tx_id"]   = mrna["attributes"].str.extract(r"ID=([^;]+)")
mrna["gene_id"] = mrna["attributes"].str.extract(r"Parent=([^;]+)")
tx2gene = mrna.set_index("tx_id")["gene_id"].to_dict()


cds = gff[gff["feature"] == "CDS"].copy()
cds["tx_id"]   = cds["attributes"].str.extract(r"Parent=([^;]+)")
cds["gene_id"] = cds["tx_id"].map(tx2gene)


protein_coding_genes = cds["gene_id"].dropna().unique()
print(f"Number of protein-coding genes: {len(protein_coding_genes)}")


per_scaffold = (
    cds.dropna(subset=["gene_id"])
       .groupby("chrom")["gene_id"]
       .nunique()
       .reset_index(name="protein_coding_gene_count")
)
print(per_scaffold)
per_scaffold.to_csv(f"protein_coding_genes_per_scaffold.tsv", 
                    sep="\t", index=False)


# Process each chromosome
for chrom, length in fai.values:
    bins = np.arange(0, length + window_size, window_size)
    bin_labels = bins[:-1]

    chr_id=chrom.split('_')[-1]
    # Filter genes on this chrom
    sub_genes = genes[genes["chrom"] == chrom]
    if sub_genes.empty:
        continue

    # Assign each gene to a bin
    sub_genes["bin"] = pd.cut(sub_genes["start"], bins=bins, labels=bin_labels, right=False)
    binned = sub_genes.groupby("bin").size().reindex(bin_labels, fill_value=0).reset_index()
    binned.columns = ["start", "gene_count"]
    binned["start"] = binned["start"].astype(int)

    # Plot
    plt.figure(figsize=(10, 3))
    plt.plot(binned["start"] / 1e6, binned["gene_count"], lw=1.2)

    plt.title(f"Gene Density – {chr_id}")
    plt.xlabel("Position (Mb)")
    plt.ylabel("Gene Count per 100kb")
    plt.ylim(0,80)

    # Add x-axis ticks every 5 Mb
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))

    plt.tight_layout()
    plt.savefig(f"Chr{chr_id}_gene_density.png", dpi=300)
    plt.close()


# Summarize gene features per chromosome
# Add gene lengths
gene_summary = genes.copy()
gene_summary["length"] = gene_summary["end"] - gene_summary["start"] + 1

# Force chromosome column to string to ensure merge compatibility
gene_summary["chrom"] = gene_summary["chrom"].astype(str)
fai["chrom"] = fai["chrom"].astype(str)

# Count gene stats by chromosome
agg = gene_summary.groupby("chrom").agg(
    gene_count=("length", "count"),
    total_gene_bp=("length", "sum"),
    mean_gene_len=("length", "mean"),
    median_gene_len=("length", "median")
).reset_index()

print(agg)
print(fai)

# Strip and ensure consistent dtype for merging
agg["chrom"] = agg["chrom"].astype(str).str.strip()
fai["chrom"] = fai["chrom"].astype(str).str.strip()

# Merge with chromosome size list (so all chromosomes show)
summary = fai.merge(agg, on="chrom", how="left").fillna(0)
summary["gene_density_per_mb"] = summary["gene_count"] / (summary["length"] / 1e6)
print(summary)
# Save
summary.to_csv(f"gene_feature_summary_by_chromosome.txt", sep="\t", index=False, float_format="%.2f")
