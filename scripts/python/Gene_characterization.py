import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict

# Set working directory
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# Inputs
gff_file = "braker.gff3"
fai_file = "primary_dedup_chr_masked_hp_sealed.fa.fai"

# Load scaffold sizes
fai = pd.read_csv(fai_file, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"])

# Load GFF3 file
gff = pd.read_csv(gff_file, sep="\t", comment="#", header=None,
                  names=["chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"])

# Extract gene, mRNA, exon, CDS
genes = gff[gff["feature"] == "gene"].copy()
mrnas = gff[gff["feature"] == "mRNA"].copy()
exons = gff[gff["feature"] == "exon"].copy()

# Helper: extract ID and Parent from attributes
def extract_id(attr):
    for item in attr.split(";"):
        if item.startswith("ID="):
            return item[3:]
    return None

def extract_parent(attr):
    for item in attr.split(";"):
        if item.startswith("Parent="):
            return item[7:]
    return None

genes["ID"] = genes["attributes"].apply(extract_id)
mrnas["ID"] = mrnas["attributes"].apply(extract_id)
mrnas["Parent"] = mrnas["attributes"].apply(extract_parent)
exons["Parent"] = exons["attributes"].apply(extract_parent)

# === Stats ===
total_genes = genes.shape[0]
total_mrnas = mrnas.shape[0]
avg_transcripts_per_gene = total_mrnas / total_genes

# Count exons per transcript
exon_counts = exons.groupby("Parent").size()
avg_exons_per_transcript = exon_counts.mean()

# Count introns per transcript (introns = exons - 1 for transcripts with ≥2 exons)
intron_counts = exon_counts[exon_counts > 1] - 1
avg_introns_per_transcript = intron_counts.mean()

# Alternative splicing: genes with >1 mRNA
transcripts_per_gene = mrnas.groupby("Parent").size()
num_alternatively_spliced_genes = (transcripts_per_gene > 1).sum()
percent_alternatively_spliced = 100 * num_alternatively_spliced_genes / total_genes

# === Summary ===
print(f"Total genes: {total_genes:,}")
print(f"Total mRNAs (transcripts): {total_mrnas:,}")
print(f"Average transcripts per gene: {avg_transcripts_per_gene:.2f}")
print(f"Average exons per transcript: {avg_exons_per_transcript:.2f}")
print(f"Average introns per transcript (with >1 exon): {avg_introns_per_transcript:.2f}")
print(f"Number of alternatively spliced genes: {num_alternatively_spliced_genes:,} ({percent_alternatively_spliced:.1f}%)")


# Count number of exons per transcript
exon_counts = exons.groupby("Parent").size()

# Get transcripts with only one exon
single_exon_transcripts = exon_counts[exon_counts == 1].index

# Subset mRNAs that are intronless
intronless_mrnas = mrnas[mrnas["ID"].isin(single_exon_transcripts)]

# Map to parent genes
intronless_gene_candidates = intronless_mrnas["Parent"]

# Option 1: Any transcript is intronless → count gene
genes_with_any_intronless_tx = set(intronless_gene_candidates)

# Option 2: All transcripts are intronless → stricter
transcripts_per_gene = mrnas.groupby("Parent")["ID"].apply(list)

strictly_intronless_genes = [
    gene for gene, tx_list in transcripts_per_gene.items()
    if set(tx_list).issubset(single_exon_transcripts)
]

# Print summary
print(f"Transcripts with only one exon: {len(single_exon_transcripts):,}")
print(f"Genes with at least one intronless transcript: {len(genes_with_any_intronless_tx):,}")
print(f"Genes with only intronless transcripts (strict): {len(strictly_intronless_genes):,}")
# print(strictly_intronless_genes)

# Histogram of exons per transcript
plt.figure(figsize=(6,4))
exon_counts.hist(bins=range(1, exon_counts.max()+2), edgecolor=None)
plt.xlabel("Exons per transcript")
plt.ylabel("Number of transcripts")
plt.xlim(0,100)
plt.title("Exon Count Distribution")
plt.savefig("Exon_count_distribution.png",dpi=600)
plt.tight_layout()


import pandas as pd
from collections import Counter

# Load the annotation file (skip comment lines)
annotations = pd.read_csv("eggnog_annotation.emapper.annotations", sep="\t", comment="#", low_memory=False)

# Total genes (assumed 1 row per gene)
total_genes = annotations.shape[0]

# Functionally annotated genes: defined by presence of any COG, GO, or KEGG info
annotated = annotations[
    (annotations["COG_category"].notna() & (annotations["COG_category"] != "-")) |
    (annotations["KEGG_ko"].notna() & (annotations["KEGG_ko"] != "-")) |
    (annotations["GOs"].notna() & (annotations["GOs"] != "-"))
]
num_annotated = annotated.shape[0]

# Genes with assigned gene names (column: "Preferred_name")
named_genes = annotations[
    (annotations["Preferred_name"].notna()) & 
    (annotations["Preferred_name"] != "-")
]
num_named = named_genes.shape[0]

# COG category breakdown
cog_counts = (
    annotations["COG_category"]
    .dropna()
    .str.cat(sep="")
    .replace("-", "")
)
cog_distribution = Counter(cog_counts)

# Top GO terms
go_series = annotations["GOs"].dropna().str.split(",")
go_counts = pd.Series([go for sublist in go_series for go in sublist]).value_counts()

# Top KEGG pathways
kegg_series = annotations["KEGG_ko"].dropna().str.split(",")
kegg_counts = pd.Series([ko for sublist in kegg_series for ko in sublist]).value_counts()

# Output
print(f"Total genes: {total_genes}")
print(f"Functionally annotated genes: {num_annotated} ({100 * num_annotated / total_genes:.2f}%)")
print(f"Genes with assigned gene names: {num_named} ({100 * num_named / total_genes:.2f}%)")

print("\nTop 10 COG categories:")
for cog, count in cog_distribution.most_common(10):
    print(f"{cog}: {count}")

print("\nTop 10 GO terms:")
print(go_counts.head(10))

print("\nTop 10 KEGG KOs:")
print(kegg_counts.head(10))






