import pandas as pd
import pybedtools
import os

# ================================
# CONFIG
# ================================
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final")

GAP_RECORDS_FILE = "gap_records_analysis.csv"
GENE_FUNCTION_BED = "genes_function.bed"
ORTHO_CATEGORIES_FILE = "all_genes_orthogroup_categories.csv"

INSERTION_BED_FILE = "cv_insertion_regions.bed"
GENES_IN_INSERTIONS_FILE = "cv_insertion_genes.csv"
ANNOTATED_OUTPUT_FILE = "cv_insertion_genes_with_orthogroups.csv"

# ================================
# 1. Load gap records and extract CV-insertions
# ================================
print("*** Loading gap records ***")
gap_df = pd.read_csv(GAP_RECORDS_FILE)

# Filter for boundary label and insertion_cv
cv_insertions = gap_df[
    (gap_df["boundary_label"] == "boundary") &
    (gap_df["event_type"] == "insertion_cv")
].copy()

print(f"Found {len(cv_insertions)} CV insertions (boundary-label filtered)")


# Create BED-style DataFrame for CV-insertions
# Create BED-style DataFrame for CV-insertions
chrom_map = {f"cf{i}": f"HiC_scaffold_{i}" for i in range(1, 100)}

cv_insertions_bed = cv_insertions[["cf_chr", "cf_pos1", "cf_pos2"]].copy()
cv_insertions_bed.columns = ["chrom", "start", "end"]
cv_insertions_bed["start"] = cv_insertions_bed["start"].astype(int)
cv_insertions_bed["end"] = cv_insertions_bed["end"].astype(int)
cv_insertions_bed["chrom"] = cv_insertions_bed["chrom"].str.strip().replace(chrom_map)

# ✅ Sort BED
cv_insertions_bed = cv_insertions_bed.sort_values(["chrom", "start"])

cv_insertions_bed.to_csv(INSERTION_BED_FILE, sep="\t", header=False, index=False)

# ================================
# 2. Intersect with gene BED
# ================================
print("*** Finding genes within CV-insertions ***")
insertions_bt = pybedtools.BedTool(INSERTION_BED_FILE)
genes_bt = pybedtools.BedTool(GENE_FUNCTION_BED)

overlaps = insertions_bt.intersect(genes_bt, wa=True, wb=True)

overlap_records = []
for f in overlaps:
    fields = f.fields
    gene_chrom = fields[3]
    gene_start = int(fields[4])
    gene_end = int(fields[5])
    gene_id = fields[6]
    gene_name = fields[7] if fields[7] else None
    strand = fields[8]
    
    overlap_records.append({
        "gene_chrom": gene_chrom,
        "gene_start": gene_start,
        "gene_end": gene_end,
        "gene_id": gene_id,
        "gene_name": gene_name,
        "strand": strand
    })

genes_in_insertions_df = pd.DataFrame(overlap_records)
print(f"Found {len(genes_in_insertions_df)} genes within CV-insertions")

genes_in_insertions_df.to_csv(GENES_IN_INSERTIONS_FILE, index=False)
print(f"*** Saved genes within insertions: {GENES_IN_INSERTIONS_FILE} ***")

# ================================
# 3. Merge with orthogroup category
# ================================
print("*** Annotating genes with orthogroup category ***")
ortho_df = pd.read_csv(ORTHO_CATEGORIES_FILE)

annotated_df = genes_in_insertions_df.merge(
    ortho_df[ortho_df["species"] == "CV"],
    on="gene_id",
    how="left"
)

# ================================
# 4. Save annotated output
# ================================
annotated_df.to_csv(ANNOTATED_OUTPUT_FILE, index=False)
print(f"*** Saved annotated genes in insertions: {ANNOTATED_OUTPUT_FILE} ***")

# Quick summary
print("\nCategory counts in annotated insertions:")
print(annotated_df["category"].value_counts(dropna=False))

print("\n*** ALL DONE ***")

