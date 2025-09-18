import pandas as pd
import numpy as np
from pathlib import Path
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")
# === USER INPUTS ===
INPUT_FILE = Path("wgs_depth_binned.txt")
fai_file = Path("primary_dedup_chr_masked_hp_sealed.fa.fai")
OUTPUT_FILE = Path("coverage.tsv")

ORIGINAL_BIN = 100_000
NEW_BIN = 200_000 

# === LOAD INPUT ===
df = pd.read_csv(INPUT_FILE, sep="\t")
print(f"Loaded {len(df)} bins.")

# === Compute new bin start ===
# Floor start to nearest new bin boundary
df["new_bin_start"] = (df["start"] // NEW_BIN) * NEW_BIN
df["new_bin_end"] = df["new_bin_start"] + NEW_BIN

# === Group by chromosome and new bins ===
agg = df.groupby(["chrom", "new_bin_start", "new_bin_end"])["depth"].mean().reset_index()

# === Rename columns to match original style ===
agg = agg.rename(columns={
    "new_bin_start": "start",
    "new_bin_end": "end"
})

# === Save to file ===
agg.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"✅ Binned coverage saved to {OUTPUT_FILE}")