import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns
import re
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster")


# Load data
df = pd.read_csv("20cv_df3_hoh_biSNP_filtered_maf05_ann_sig5_sam8_functional.txt", sep="\t")

# Extract sample columns
sample_cols = df.columns[5:]  # assuming first 5 are metadata

# Define L and S groups
L_samples = [col for col in sample_cols if col.endswith("L_st_dr_filtered")]
S_samples = [col for col in sample_cols if col.endswith("S_st_dr_filtered")]

# Normalize genotype formats
def parse_genotype(raw):
    gt = str(raw).split(";")[0]
    if gt in ["0/0", "0|0"]: return "0/0"
    elif gt in ["0/1", "1/0", "0|1", "1|0"]: return "0/1"
    elif gt in ["1/1", "1|1"]: return "1/1"
    else: return "NA"

# Count genotypes
results = []
for idx, row in df.iterrows():
    for group_name, samples in [("L", L_samples), ("S", S_samples)]:
        genotypes = [parse_genotype(row[s]) for s in samples]
        counts = Counter(genotypes)
        results.append({
            "Locus": f"{row['#Chromosome']}:{row['position']}",
            "Group": group_name,
            "0/0": counts["0/0"],
            "0/1": counts["0/1"],
            "1/1": counts["1/1"],
            "Missing": counts["NA"]
        })

# Convert to DataFrame and save
result_df = pd.DataFrame(results)
result_df.to_csv("genotype_counts_by_group.tsv", sep="\t", index=False)
print(result_df)