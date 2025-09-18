import pandas as pd
import numpy as np
from collections import defaultdict
import os

# === SETTINGS ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

INPUT_FILE = "coords_large_sv.tsv"
FAI_FILE = "primary_dedup_chr_masked_hp_sealed.fa.fai"
OUTPUT_FILE = "SD_density_by_type_100kb.tsv"

WINDOW_SIZE = 500_000
MIN_LENGTH = 10_000

# === Load chromosome sizes from .fai ===
print("Loading .fai file...")
fai_df = pd.read_csv(FAI_FILE, sep="\t", header=None,usecols=[0, 1], names=["chr", "length"])
chrom_sizes = dict(zip(fai_df["chr"], fai_df["length"]))
print(f"Loaded {len(chrom_sizes)} chromosomes.")

# === Make fixed-size windows ===
print(f"Generating {WINDOW_SIZE}-bp windows...")
windows = []
for chrom, size in chrom_sizes.items():
    for start in range(0, size, WINDOW_SIZE):
        end = min(start + WINDOW_SIZE, size)
        windows.append((chrom, start, end))

windows_df = pd.DataFrame(windows, columns=["chr", "start", "end"])
windows_df["Interchromosomal_bp"] = 0
windows_df["SD_bp"] = 0
windows_df["TD_bp"] = 0
print(f"Total windows created: {len(windows_df)}.")

# === Load SV data ===
print("Reading alignment file...")
sv_df = pd.read_csv(INPUT_FILE, sep="\t")
print(f"Total alignments in file: {len(sv_df)}.")

# Filter criteria
sv_df = sv_df[(sv_df["sv_candidate"] == True) & (sv_df["alignment_length"] >= MIN_LENGTH)]
print(f"Retained after filtering: {len(sv_df)} alignments.")

# === Category mapping ===
CATEGORY_MAP = {
    "Segmental Duplication (Inter-chromosomal)": "Interchromosomal_bp",
    "Segmental Duplication": "SD_bp",
    "Tandem Duplication": "TD_bp"
}

# === Collect intervals (both REF and QRY sides) ===
all_intervals = []
for _, row in sv_df.iterrows():
    category = row["category"].strip()
    if category not in CATEGORY_MAP:
        continue
    out_col = CATEGORY_MAP[category]

    # Add REF side
    all_intervals.append((row["ref_chr"], row["ref_start"], row["ref_end"], out_col))
    # Add QRY side
    all_intervals.append((row["qry_chr"], row["qry_start"], row["qry_end"], out_col))

intervals_df = pd.DataFrame(all_intervals, columns=["chr", "start", "end", "type"])
intervals_df = intervals_df[intervals_df["chr"].isin(chrom_sizes)]
print(f"Total intervals (REF+QRY sides): {len(intervals_df)}.")

# === Bin intervals into windows ===
print("Binning intervals into windows...")
windows_by_chr = {chrom: df for chrom, df in windows_df.groupby("chr")}

for chrom, group in intervals_df.groupby("chr"):
    if chrom not in windows_by_chr:
        continue
    chrom_windows = windows_by_chr[chrom]

    for _, interval in group.iterrows():
        iv_start, iv_end, iv_type = interval["start"], interval["end"], interval["type"]

        overlaps = (chrom_windows["end"] > iv_start) & (chrom_windows["start"] < iv_end)
        overlap_windows = chrom_windows.loc[overlaps]

        for idx, win in overlap_windows.iterrows():
            overlap_start = max(iv_start, win["start"])
            overlap_end = min(iv_end, win["end"])
            overlap_len = overlap_end - overlap_start
            if overlap_len > 0:
                windows_df.loc[idx, iv_type] += overlap_len

print("Finished binning.")

# === Save final TSV ===
windows_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"✅ Binned SD density track written to: {OUTPUT_FILE}")