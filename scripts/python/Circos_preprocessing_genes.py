import pandas as pd
import numpy as np
from pathlib import Path
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")
# === USER INPUTS ===
gene_bed = Path("genes_function.bed")
fai_file = Path("primary_dedup_chr_masked_hp_sealed.fa.fai")
window_size = 100_000
output_tsv = Path("gene_density.tsv")


# === STEP 1: Load gene BED ===
genes = pd.read_csv(gene_bed, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"])
genes["chrom"] = genes["chrom"].astype(str).str.strip()
genes["start"] = genes["start"].astype(int)
genes["end"]   = genes["end"].astype(int)
print(genes)
# === STEP 2: Load genome sizes and define bins ===
chrom_sizes = pd.read_csv(fai_file, sep="\t", header=None, usecols=[0, 1],
                          names=["chrom", "size"])
chrom_sizes["chrom"] = chrom_sizes["chrom"].astype(str).str.strip()

bins = []
for _, row in chrom_sizes.iterrows():
    chrom, size = row["chrom"], row["size"]
    edges = np.arange(0, size + window_size, window_size)
    bins += [
        {"chrom": chrom,
         "bin_start": int(edges[i]),
         "bin_end":   int(edges[i+1])}
        for i in range(len(edges) - 1)
    ]
bin_df = pd.DataFrame(bins)

# === STEP 3: Assign each gene to overlapping bins ===
def find_bins(g):
    sel = bin_df[
        (bin_df["chrom"] == g["chrom"]) &
        (g["start"] < bin_df["bin_end"]) &
        (g["end"] > bin_df["bin_start"])
    ]
    return sel.index.tolist()

gene_to_bins = genes.apply(find_bins, axis=1)
valid_mask = gene_to_bins.str.len() > 0
genes_valid = genes[valid_mask].reset_index(drop=True)
gene_to_bins = gene_to_bins[valid_mask].reset_index(drop=True)

# === STEP 4: Expand gene-bin relationships and count genes per bin ===
genes_exp = genes_valid.loc[genes_valid.index.repeat(gene_to_bins.str.len())].copy()
genes_exp["bin_index"] = [i for sublist in gene_to_bins for i in sublist]
bin_df = bin_df.reset_index().rename(columns={"index": "bin_index"})

# merge with bin info (preserve chrom from bin)
genes_exp = genes_exp.merge(bin_df, on="bin_index", suffixes=("", "_bin"))
genes_exp["chrom"] = genes_exp["chrom_bin"]

# === STEP 5: Count genes per bin
density = (
    genes_exp.groupby(["chrom", "bin_start", "bin_end"])
    .size()
    .reset_index(name="gene_count")
)

# === STEP 6: Save output
density.to_csv(output_tsv, sep="\t", index=False)