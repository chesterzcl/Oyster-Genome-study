import re
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/comparative_genomics"
collinearity_file = os.path.join(dir, "CV_Mgiga.collinearity")
functional_bed_file = os.path.join(dir, "genes_function.bed")  # your BED-like file

alignment_blocks = []
current_block = None

with open(collinearity_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith("## Alignment"):
            if current_block:
                alignment_blocks.append(current_block)

            block_id = int(re.search(r"Alignment\s+(\d+):", line).group(1))
            score = float(re.search(r"score=([\d\.]+)", line).group(1))
            num_pairs = int(re.search(r"N=(\d+)", line).group(1))

            tokens = re.split(r'\s+', line)
            scaffold_info = tokens[-2:]

            current_block = {
                "block_id": block_id,
                "score": score,
                "num_pairs": num_pairs,
                "genes": [],
                "scaffolds": [s.strip() for s in scaffold_info]
            }

        elif not line.startswith("#"):
            parts = line.split()
            if len(parts) >= 5:
                gene1 = parts[2].strip()
                gene2 = parts[3].strip()
                current_block["genes"].append({
                    "block_id": current_block["block_id"],
                    "gene1": gene1,
                    "gene2": gene2,
                    "evalue": 0.0  # MCScanX may not output e-value, placeholder
                })

if current_block:
    alignment_blocks.append(current_block)

# === Create DataFrames ===
block_summary = pd.DataFrame([{
    "block_id": b["block_id"],
    "score": b["score"],
    "num_pairs": b["num_pairs"],
    "scaffold1": b["scaffolds"][0],
    "scaffold2": b["scaffolds"][1],
} for b in alignment_blocks])

gene_pairs = pd.DataFrame([
    gene for block in alignment_blocks for gene in block["genes"]
])

# Save summary files
block_summary.to_csv(os.path.join(dir, "mcscanx_block_summary.csv"), index=False)
gene_pairs.to_csv(os.path.join(dir, "mcscanx_gene_pairs.csv"), index=False)

# === Load Functional BED file ===
columns = ["chrom", "start", "end", "gene_id", "symbol", "biotype"]
gene_bed = pd.read_csv(functional_bed_file, sep="\t", names=columns)

# === Merge to get coordinates for C. virginica genes ===
merged = gene_pairs.merge(gene_bed, left_on="gene1", right_on="gene_id")

# === Group by block to get BED coordinates ===
block_bed = merged.groupby("block_id").agg(
    chrom=("chrom", lambda x: x.mode().iloc[0]),  # dominant scaffold
    start=("start", "min"),
    end=("end", "max")
).reset_index()

block_bed["name"] = "synteny_block_" + block_bed["block_id"].astype(str)

# === Output BED ===
bed_output = block_bed[["chrom", "start", "end", "name"]]


def extract_chr_num(chrom):
    match = re.search(r'HiC_scaffold_(\d+)', chrom)
    return int(match.group(1)) if match else float('inf')

bed_output["chr_num"] = bed_output["chrom"].apply(extract_chr_num)
bed_output = bed_output.sort_values(by=["chr_num", "start"])
bed_output = bed_output.drop(columns=["chr_num"])
print(bed_output)

# Calculate total length of all synteny blocks
bed_output["length"] = bed_output["end"] - bed_output["start"]
total_length = bed_output["length"].sum()
print(f"Total synteny block length: {total_length:,} bp")

bed_output.to_csv(os.path.join(dir, "synteny_blocks_cv.bed"), sep="\t", header=False, index=False)


syntenic_genes = set(gene_pairs["gene1"].unique())
gene_bed["is_syntenic"] = gene_bed["gene_id"].isin(syntenic_genes)
gene_bed["chr_num"] = gene_bed["chrom"].apply(extract_chr_num)
gene_bed = gene_bed.sort_values(by=["chr_num", "start"])
gene_bed = gene_bed.drop(columns=["chr_num"])

syntenic_df = gene_bed[gene_bed["is_syntenic"]].copy().drop(columns=["is_syntenic"])
nonsyntenic_df = gene_bed[~gene_bed["is_syntenic"]].copy().drop(columns=["is_syntenic"])
# print(syntenic_df)
# print(nonsyntenic_df)
# print(gene_bed)

syntenic_df.to_csv(os.path.join(dir,"syntenic_genes.bed"), sep="\t", index=False, header=False)
nonsyntenic_df.to_csv(os.path.join(dir,"non_syntenic_genes.bed"), sep="\t", index=False, header=False)


from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

# Use mapped orthologs if needed
go_syntenic = gp.profile(organism='hsapiens', query=syntenic_df["symbol"].dropna().tolist())
go_nonsyntenic = gp.profile(organism='hsapiens', query=nonsyntenic_df["symbol"].dropna().tolist())

# Filter for Biological Process
go_syntenic_bp = go_syntenic[go_syntenic['source'] == 'GO:BP']
go_nonsyntenic_bp = go_nonsyntenic[go_nonsyntenic['source'] == 'GO:BP']

# Save
go_syntenic_bp.to_csv(os.path.join(dir,"go_enrichment_syntenic.csv"), index=False)
go_nonsyntenic_bp.to_csv(os.path.join(dir,"go_enrichment_nonsyntenic.csv"), index=False)


def plot_go(go_df, title):
    top = go_df.sort_values("p_value").head(10)
    plt.figure(figsize=(8, 6))
    plt.barh(top["name"], -np.log10(top["p_value"]))
    plt.xlabel("-log10(p-value)")
    plt.title(title)
    plt.tight_layout()
    # plt.show()

# plot_go(go_syntenic_bp, "Top GO:BP Enriched in Syntenic Genes")
# plot_go(go_nonsyntenic_bp, "Top GO:BP Enriched in Non-Syntenic Genes")

# Load high- and low-variance ATAC peak files
high_atac = pd.read_csv(os.path.join(dir,"ATAC_high_variance_peaks.bed"), sep="\t", header=None, names=["chrom", "start", "end"])
low_atac = pd.read_csv(os.path.join(dir,"ATAC_low_variance_peaks.bed"), sep="\t", header=None, names=["chrom", "start", "end"])

# Function to find overlaps between two BED-like DataFrames
def find_overlaps(atac_df, synteny_df):
    overlaps = []
    for _, atac_row in atac_df.iterrows():
        matches = synteny_df[
            (synteny_df["chrom"] == atac_row["chrom"]) &
            (synteny_df["start"] < atac_row["end"]) &
            (synteny_df["end"] > atac_row["start"])
        ]
        for _, match in matches.iterrows():
            overlaps.append({
                "peak_chrom": atac_row["chrom"],
                "peak_start": atac_row["start"],
                "peak_end": atac_row["end"],
                "block_id": match["name"],
                "block_start": match["start"],
                "block_end": match["end"]
            })
    return pd.DataFrame(overlaps)

# Run overlap detection
high_in_synteny = find_overlaps(high_atac, bed_output)
low_in_synteny = find_overlaps(low_atac, bed_output)

# Save results
high_in_synteny.to_csv(os.path.join(dir,"ATAC_high_in_synteny.csv"), index=False)
low_in_synteny.to_csv(os.path.join(dir,"ATAC_low_in_synteny.csv"), index=False)

# Count summary
print(f"High-variance peaks in synteny: {len(high_in_synteny)} / {len(high_atac)}")
print(f"Low-variance peaks in synteny: {len(low_in_synteny)} / {len(low_atac)}")

from scipy.stats import fisher_exact

# Build the table
table = [[10172, 3263],
         [64343, 12468]]

# Run Fisher's exact test
odds_ratio, p_value = fisher_exact(table)

print(f"Odds ratio: {odds_ratio:.3f}")
print(f"P-value: {p_value:.3e}")

