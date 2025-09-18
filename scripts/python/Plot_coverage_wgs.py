import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === Paths ===
hifi_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/hifi_depth_binned.txt"
illumina_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/wgs_depth_binned.txt"
output_dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/coverage_plot"

# === Load and normalize ===
hifi = pd.read_csv(hifi_file, sep="\t")
illumina = pd.read_csv(illumina_file, sep="\t")

# Rename for clarity
hifi = hifi.rename(columns={"depth": "hifi_depth"})
illumina = illumina.rename(columns={"depth": "illumina_depth"})

# Merge on chrom/start
merged = pd.merge(hifi, illumina, on=["chrom", "start"])

# Normalize coverage
merged["hifi_norm"] = merged["hifi_depth"] / merged["hifi_depth"].mean()
merged["illumina_norm"] = merged["illumina_depth"] / merged["illumina_depth"].mean()

# === Step: Define half-coverage threshold ===
threshold = 0.65  # adjust as needed

# Extract bins with both HiFi and Illumina normalized coverage below threshold
half_coverage = merged[
    (merged["hifi_norm"] < threshold) &
    (merged["illumina_norm"] < threshold)
].copy()

# Optionally add end coordinate (if bin size is known, e.g., 10kb)
bin_size = merged["start"].diff().median()  # or hardcode if known
half_coverage["end"] = half_coverage["start"] + bin_size

# Save to BED-like format
half_coverage_bed = half_coverage[["chrom", "start", "end"]]
half_coverage_bed.to_csv(f"{output_dir}/half_coverage_regions.bed", sep="\t", index=False, header=False)

print(f"Extracted {len(half_coverage)} half-coverage bins to: half_coverage_regions.bed")
# === Plot per chromosome ===
for chrom in merged['chrom'].unique():
    sub = merged[merged['chrom'] == chrom].copy()
    chr_id = chrom.split('_')[-1]

    # Convert to Mb
    sub["start_MB"] = sub["start"] / 1e6
    max_mb = sub["start_MB"].max()

    plt.figure(figsize=(12, 4))
    plt.plot(sub["start_MB"], sub["hifi_norm"], label="HiFi (normalized)", lw=1.5)
    plt.plot(sub["start_MB"], sub["illumina_norm"], label="Illumina (normalized)", lw=1.5)

    plt.title(f"Normalized Coverage on Chr{chr_id}")
    plt.xlabel("Position (MB)")
    plt.ylabel("Normalized Coverage")
    plt.ylim(0, 4)  # adjust as needed
    plt.xticks(np.arange(0, max_mb + 5, 5))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/Chr{chr_id}_normalized_coverage.png", dpi=300)
    plt.close()

print("Normalized coverage plots saved in:", output_dir)