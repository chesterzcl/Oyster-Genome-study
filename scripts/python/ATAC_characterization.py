import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Set working directory ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")
#########################################################
# SETTINGS
#########################################################

PROMOTER_UPSTREAM = 1000
PROMOTER_DOWNSTREAM = 100

INPUT_PEAK_FILE = "consensus_counts_matrix.txt"
INPUT_GFF_FILE = "braker_with_name_and_description.gff3"
OUTPUT_FILE = "ATAC_peaks_annotated_with_labels.txt"

#########################################################
# 1️⃣ Load ATAC peaks
#########################################################

print("Loading peaks...")
peaks = pd.read_csv(INPUT_PEAK_FILE, sep="\t", usecols=["chr", "start", "end"])
peaks["peak_id"] = peaks.index
print(f"Loaded {len(peaks)} peaks.")

#########################################################
# 2️⃣ Load gene annotations from GFF3
#########################################################

print("Loading GFF3 genes...")
gff = pd.read_csv(
    INPUT_GFF_FILE,
    sep="\t",
    header=None,
    comment="#",
    names=[
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ]
)

# Filter for gene entries
genes = gff[gff["type"] == "gene"].copy()
genes["gene_id"] = genes["attributes"].str.extract(r'ID=([^;]+)')
genes["gene_name"] = genes["attributes"].str.extract(r'Name=([^;]+)')
genes["gene_name"] = genes["gene_name"].fillna(genes["gene_id"])
genes = genes[["seqid", "start", "end", "strand", "gene_id", "gene_name"]]
print(f"Loaded {len(genes)} genes.")

#########################################################
# 3️⃣ Define promoter regions
#########################################################

print("Defining promoter regions...")
def define_promoter(row, upstream, downstream):
    if row["strand"] == "+":
        tss = row["start"]
        prom_start = max(0, tss - upstream)
        prom_end = tss + downstream
    else:
        tss = row["end"]
        prom_start = max(0, tss - downstream)
        prom_end = tss + upstream
    return pd.Series([prom_start, prom_end])

genes[["promoter_start", "promoter_end"]] = genes.apply(
    define_promoter,
    axis=1,
    upstream=PROMOTER_UPSTREAM,
    downstream=PROMOTER_DOWNSTREAM
)

#########################################################
# 4️⃣ Classify peaks: promoter / gene_body / distal
#########################################################

print("Classifying promoter overlaps...")
# Merge for promoter overlap
merged_promoter = peaks.merge(
    genes,
    left_on="chr",
    right_on="seqid",
    how="left"
)

promoter_overlap = merged_promoter[
    (merged_promoter["start_x"] < merged_promoter["promoter_end"]) &
    (merged_promoter["end_x"] > merged_promoter["promoter_start"])
]

promoter_peak_ids = promoter_overlap["peak_id"].unique()
promoter_gene_map = promoter_overlap.groupby("peak_id")["gene_id"].first()

print("Classifying gene-body overlaps...")
# Merge for gene-body overlap
merged_gene = peaks.merge(
    genes,
    left_on="chr",
    right_on="seqid",
    how="left"
)

gene_body_overlap = merged_gene[
    (merged_gene["start_x"] < merged_gene["end_y"]) &
    (merged_gene["end_x"] > merged_gene["start_y"])
]

gene_body_peak_ids = gene_body_overlap["peak_id"].unique()
gene_body_gene_map = gene_body_overlap.groupby("peak_id")["gene_id"].first()

# Initialize all as distal
peaks["location"] = "distal"
peaks["label"] = None

# Assign gene_body
peaks.loc[peaks["peak_id"].isin(gene_body_peak_ids), "location"] = "gene_body"
peaks.loc[peaks["peak_id"].isin(gene_body_peak_ids), "label"] = peaks["peak_id"].map(gene_body_gene_map)

# Assign promoter (priority)
peaks.loc[peaks["peak_id"].isin(promoter_peak_ids), "location"] = "promoter"
peaks.loc[peaks["peak_id"].isin(promoter_peak_ids), "label"] = peaks["peak_id"].map(promoter_gene_map)

#########################################################
# 5️⃣ For distal peaks: assign upstream-downstream gene IDs
#########################################################

print("Labeling distal peaks...")
distal_peaks = peaks[peaks["location"] == "distal"].copy()

distal_labels = []

for idx, row in distal_peaks.iterrows():
    chrom = row["chr"]
    p_start = row["start"]
    p_end = row["end"]
    
    genes_chr = genes[genes["seqid"] == chrom]
    
    upstream_genes = genes_chr[genes_chr["end"] < p_start]
    downstream_genes = genes_chr[genes_chr["start"] > p_end]
    
    if not upstream_genes.empty:
        upstream_gene = upstream_genes.iloc[upstream_genes["end"].argmax()]["gene_id"]
    else:
        upstream_gene = "NA"
        
    if not downstream_genes.empty:
        downstream_gene = downstream_genes.iloc[downstream_genes["start"].argmin()]["gene_id"]
    else:
        downstream_gene = "NA"
    
    intergenic_label = f"{upstream_gene}-{downstream_gene}"
    distal_labels.append(intergenic_label)

peaks.loc[peaks["location"] == "distal", "label"] = distal_labels

#########################################################
# 6️⃣ Save final annotated peaks
#########################################################

print(f"Saving to {OUTPUT_FILE}...")
peaks[["chr", "start", "end", "location", "label"]].to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

print("Done!")






