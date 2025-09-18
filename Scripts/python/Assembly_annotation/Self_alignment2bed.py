import pandas as pd
import os


os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/residual_haplotig")

# Load classified coords file

df = pd.read_csv("coords_large_sv.tsv", sep="\t")

# Set minimum alignment size
MIN_LEN = 10000

# Filter: keep rows where both ref and qry regions are long enough
filtered_df = df[
    ((df["ref_end"] - df["ref_start"]) >= MIN_LEN) &
    ((df["qry_end"] - df["qry_start"]) >= MIN_LEN)
].copy()

# Extract ref and qry regions
ref_bed = filtered_df[["ref_chr", "ref_start", "ref_end"]].copy()
qry_bed = filtered_df[["qry_chr", "qry_start", "qry_end"]].copy()

# Rename columns to standard BED format
ref_bed.columns = ["chrom", "start", "end"]
qry_bed.columns = ["chrom", "start", "end"]

# Convert to 0-based BED (BED format: start is 0-based, end is 1-based)
ref_bed["start"] = ref_bed["start"] - 1
qry_bed["start"] = qry_bed["start"] - 1

# Combine, remove duplicates, and sort
bed_df = pd.concat([ref_bed, qry_bed]).drop_duplicates()
bed_df = bed_df.sort_values(by=["chrom", "start", "end"])

# Save output BED file
bed_df.to_csv("regions_large.bed", sep="\t", header=False, index=False)

