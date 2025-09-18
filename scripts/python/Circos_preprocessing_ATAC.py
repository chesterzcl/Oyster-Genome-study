import pandas as pd
import numpy as np
from pathlib import Path
import os 

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === USER INPUTS ===
atac_bed = Path("consensus_counts_matrix.txt")      # BED format: chrom, start, end
fai_file = Path("primary_dedup_chr_masked_hp_sealed.fa.fai")
window_size = 200_000
output_tsv = Path("atac_peak_density_200k.tsv")

# === STEP 1: Load ATAC BED ===
atac = pd.read_csv(atac_bed,sep="\t",header=0,usecols=[0, 1, 2],
                   names=["chrom","start","end"])

atac["chrom"] = atac["chrom"].astype(str).str.strip()
atac["start"] = atac["start"].astype(int)
atac["end"]   = atac["end"].astype(int)
atac["length"] = atac["end"] - atac["start"]

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

# === STEP 3: Assign peaks to overlapping bins ===
def find_bins(p):
    sel = bin_df[
        (bin_df["chrom"] == p["chrom"]) &
        (p["start"] < bin_df["bin_end"]) &
        (p["end"] > bin_df["bin_start"])
    ]
    return sel.index.tolist()

peak_to_bins = atac.apply(find_bins, axis=1)
valid_mask = peak_to_bins.str.len() > 0
atac_valid = atac[valid_mask].reset_index(drop=True)
peak_to_bins = peak_to_bins[valid_mask].reset_index(drop=True)

# === STEP 4: Expand peak-bin relationships and compute overlaps ===
atac_exp = atac_valid.loc[atac_valid.index.repeat(peak_to_bins.str.len())].copy()
atac_exp["bin_index"] = [i for sublist in peak_to_bins for i in sublist]
bin_df = bin_df.reset_index().rename(columns={"index": "bin_index"})
atac_exp = atac_exp.merge(bin_df, on="bin_index", suffixes=("", "_bin"))
atac_exp["chrom"] = atac_exp["chrom_bin"]

def compute_overlap(row):
    return max(0, min(row["end"], row["bin_end"]) - max(row["start"], row["bin_start"]))
atac_exp["overlap"] = atac_exp.apply(compute_overlap, axis=1)

# === STEP 5: Sum peak coverage per bin
summary = (
    atac_exp.groupby(["chrom", "bin_start", "bin_end"])["overlap"]
    .sum()
    .reset_index()
    .rename(columns={"overlap": "atac_bp"})
)
summary["atac_pct"] = summary["atac_bp"] / window_size * 100

# === STEP 6: Save output
summary.to_csv(output_tsv, sep="\t", index=False)
print(f"ATAC peak density saved to: {output_tsv.resolve()}")