import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

bed_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/repeat_analysis/all_repeats.bed"
fai_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/primary_dedup_chr.fa.fai"
out_dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/repeat_analysis/"
window_size = 500000
focus_classes = ["DNA", "LTR", "LINE", "SINE","RC"]

legend_labels = {
    "DNA": "DNA Transposons",
    "LTR": "LTR Retrotransposons",
    "LINE": "LINE Elements",
    "SINE": "SINE Elements",
    "RC": "Rolling Circle (Helitron)"
}

df = pd.read_csv(bed_file, sep="\t", names=["chrom", "start", "end", "repeat_id", "dot", "class/family"])
df["chrom"] = df["chrom"].astype(str).str.strip()
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)
df["class"] = df["class/family"].str.extract(r"^([^/]+)")
df = df[df["class"].isin(focus_classes)]
df["length"] = df["end"] - df["start"]

# === Step 2: Load chromosome sizes and create bins ===
chrom_sizes = pd.read_csv(fai_file, sep="\t", header=None, usecols=[0, 1], names=["chrom", "size"])
chrom_sizes["chrom"] = chrom_sizes["chrom"].astype(str).str.strip()

bin_records = []
for _, row in chrom_sizes.iterrows():
    chrom = row["chrom"]
    size = row["size"]
    bins = np.arange(0, size, window_size)
    bins = np.append(bins, size)  # ensure the last bin reaches the end
    for i in range(len(bins) - 1):
        bin_records.append({
            "chrom": chrom,
            "bin_start": bins[i],
            "bin_end": bins[i + 1],
            "bin_mid": (bins[i] + bins[i + 1]) // 2
        })

bin_df = pd.DataFrame(bin_records)

# === Step 3: Assign each repeat to overlapping bins ===
def assign_bins(row):
    overlaps = bin_df[
        (bin_df["chrom"] == row["chrom"]) &
        (row["start"] < bin_df["bin_end"]) &
        (row["end"] > bin_df["bin_start"])
    ]
    return overlaps.index.tolist()

df = df.reset_index(drop=True)
repeat_to_bins = df.apply(assign_bins, axis=1)

# Filter out repeats with no bin match
valid_mask = repeat_to_bins.str.len() > 0
df_valid = df[valid_mask].copy()
repeat_to_bins = repeat_to_bins[valid_mask].reset_index(drop=True)

# Expand repeat-bin relationships
repeat_expanded = df_valid.loc[df_valid.index.repeat(repeat_to_bins.str.len())].copy()
bin_indices = [i for sublist in repeat_to_bins for i in sublist]
repeat_expanded["bin_index"] = bin_indices

# Join with bin info
bin_info = bin_df.loc[repeat_expanded["bin_index"]].reset_index(drop=True).drop(columns=["chrom"])
repeat_expanded = repeat_expanded.join(bin_info)

# === Step 4: Compute observed density per bin per class ===
observed = (
    repeat_expanded
    .groupby(["chrom", "bin_start", "bin_end", "bin_mid", "class"])["length"]
    .sum()
    .reset_index()
)
observed["density"] = observed["length"] / window_size * 100

# === Step 5: Build complete bin × class matrix and merge ===
class_df = pd.DataFrame({'class': focus_classes})
complete_bins = bin_df.assign(key=1).merge(class_df.assign(key=1), on="key").drop(columns=["key"])
grouped = pd.merge(
    complete_bins,
    observed,
    on=["chrom", "bin_start", "bin_end", "bin_mid", "class"],
    how="left"
)
grouped["density"] = grouped["density"].fillna(0)
grouped["class_pretty"] = grouped["class"].map(legend_labels)


# === Step 6: Plot per chromosome ===
for chrom in chrom_sizes["chrom"]:
    chr_data = grouped[grouped["chrom"] == chrom]
    if chr_data.empty:
        print(f"Skipping {chrom} — no bins.")
        continue
    chr_id = chrom.split('_')[-1]
    plt.figure(figsize=(12, 5))
    sns.lineplot(
        x=chr_data["bin_mid"]/1000000,
        y="density",
        hue="class_pretty",
        data=chr_data,
        linewidth=1.2,
        legend="full"
    )
    plt.title(f"Transposable element class density-Chr{chr_id}")
    plt.xlabel("Genomic Position (Mb)")
    plt.ylabel("% of Window Masked")
    plt.ylim(0,40)
    max_x = chr_data["bin_mid"].max()/1_000_000
    plt.xticks(np.arange(0,max_x+5,5))
    plt.legend(title="Class",bbox_to_anchor=(1.01, 1),loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,f"TE_class_density_Chr{chr_id}_{window_size}.png"),dpi=300)
    plt.close()