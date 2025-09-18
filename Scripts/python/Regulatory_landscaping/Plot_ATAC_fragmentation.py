import os
import matplotlib.pyplot as plt
import seaborn as sns

# --------------------------------------------------------
OUTPUT_FILE = "fragment_length_distribution_smoothed.png"
MAX_FRAGMENT_LENGTH = 500
DPI = 600
# --------------------------------------------------------

# === SETTINGS ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/fragment_lengths")

# === FIND INPUT FILES ===
files = [f for f in os.listdir() if f.endswith(".txt")]

if not files:
    print("No .txt files found in this directory. Exiting.")
    exit(1)

# === Initialize Plot ===
plt.figure(figsize=(10, 6))

for file in files:
    sample_name = file.replace("_frag_lengths.txt", "")
    # Clean to get just the first part
    short_name = sample_name.split("_")[0]
    print(f"Processing {short_name}")

    # Load fragment lengths
    with open(file) as f:
        lengths = [int(line.strip()) for line in f if line.strip().isdigit()]
    lengths = [l for l in lengths if 0 < l < MAX_FRAGMENT_LENGTH]

    if not lengths:
        print(f"[WARNING] No valid fragments in {file} after filtering. Skipping.")
        continue

    # Smoothed density plot with cleaned label
    sns.kdeplot(
        lengths,
        bw_adjust=1.0,
        label=short_name,
        alpha=0.6,
        fill=False,
        common_norm=False,
    )

# === Finalize plot ===
plt.xlabel("Fragment Length (bp)")
plt.ylabel("Density")
# plt.title("ATAC-seq Fragment Size Distribution (Smoothed KDE, All Samples)")
plt.xlim(0, MAX_FRAGMENT_LENGTH)
plt.legend()
plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=DPI)
plt.close()

print(f"\nCombined smoothed plot saved as {OUTPUT_FILE}")