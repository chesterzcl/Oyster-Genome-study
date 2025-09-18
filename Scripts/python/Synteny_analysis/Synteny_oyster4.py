import pandas as pd
import pyranges as pr
from glob import glob
from collections import defaultdict
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/oyster4/")

synteny_files = {
    "oe": "synteny_block_summary_oe.csv",
    "mg": "synteny_block_summary_mg.csv",
    "ma": "synteny_block_summary_ma.csv"
}

synteny_dict = {}
all_breakpoints = defaultdict(set)

for species, filepath in synteny_files.items():
    df = pd.read_csv(filepath, sep=",")
    df = df.rename(columns={
        "chrA": "Chromosome", "startA": "Start", "endA": "End"
    })
    df["Species"] = species
    synteny_dict[species] = df

    # Collect breakpoints for later partitioning
    for _, row in df.iterrows():
        all_breakpoints[row["Chromosome"]].update([row["Start"], row["End"]])

# Step 2: Create minimal non-overlapping segments per chromosome
segment_list = []
for chrom, breaks in all_breakpoints.items():
    sorted_breaks = sorted(breaks)
    for i in range(len(sorted_breaks) - 1):
        segment_list.append({
            "Chromosome": chrom,
            "Start": sorted_breaks[i],
            "End": sorted_breaks[i+1]
        })

segments_df = pd.DataFrame(segment_list)
segments = pr.PyRanges(segments_df)

# Step 3: Annotate presence/orientation for each species
def get_presence(segment_row, species_df):
    seg_chr = segment_row["Chromosome"]
    seg_start = segment_row["Start"]
    seg_end = segment_row["End"]

    # Get overlapping blocks
    overlaps = species_df[
        (species_df["Chromosome"] == seg_chr) &
        (species_df["Start"] < seg_end) &
        (species_df["End"] > seg_start)
    ]

    if overlaps.empty:
        return "non-syntenic"
    else:
        orientations = set(overlaps["orientation"])
        if "minus" in orientations and "plus" in orientations:
            return "inverted"  # or "mixed", depending on your biological interpretation
        elif "minus" in orientations:
            return "inverted"
        else:
            return "syntenic"

# Step 4: Apply to all intervals
annotated_segments = []
for _, seg in segments.df.iterrows():
    row = {
        "Chromosome": seg["Chromosome"],
        "Start": seg["Start"],
        "End": seg["End"]
    }
    for sp, df in synteny_dict.items():
        row[f"{sp}_status"] = get_presence(seg, df)
    annotated_segments.append(row)

result_df = pd.DataFrame(annotated_segments)

# Step 5: Classify segment by shared pattern
def classify(row):
    oe = row["oe_status"]
    mg = row["mg_status"]
    ma = row["ma_status"]

    # Collapse magallana status
    if mg == ma:
        maga = mg
    else:
        return "complex_or_conflict"

    # CASE 1: OE syntenic
    if oe == "syntenic":
        if maga == "syntenic":
            return "conserved"
        elif maga == "inverted":
            return "magallana_inversion"
        elif maga == "non-syntenic":
            return "magallana_innovation"

    # CASE 2: OE non-syntenic
    elif oe == "non-syntenic":
        if maga == "syntenic":
            return "Crassostreinae_innovation"
        elif maga == "inverted":
            return "Crassostreinae_innovation"  # Optional: you can merge or keep
        elif maga == "non-syntenic":
            return "cv_innovation"

    # CASE 3: OE inverted
    elif oe == "inverted":
        if maga == "syntenic":
            return "oe_inversion"
        elif maga == "inverted":
            return "cv_inversion"
        elif maga == "non-syntenic":
            return "magallana_innovation"

    return "complex_or_conflict"


result_df["classification"] = result_df.apply(classify, axis=1)

# Frequency summary
freq_summary = result_df["classification"].value_counts()
print(freq_summary)

# Step 6: Output
result_df.to_csv("cv_synteny_segment_evolution_summary.tsv", sep="\t", index=False)
print("Saved: cv_synteny_segment_evolution_summary.tsv")


import pandas as pd
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
import seaborn as sns
import numpy as np

# === Load your merged block length TSV ===
file_path = "synteny_block_lengths.tsv"  # Replace with actual path
df = pd.read_csv(file_path, sep="\t")

# === Reshape for violin plot ===
df_long = pd.melt(
    df,
    id_vars=["comparison"],
    value_vars=["block_length_spe1", "block_length_cv"],
    var_name="species",
    value_name="block_length"
)

# Rename species for plot labels
df_long["species"] = df_long["species"].replace({
    "block_length_spe1": "Other oysters",
    "block_length_cv": "C. virginica"
})

# Log-transform
df_long["log10_length"] = np.log10(df_long["block_length"] + 1)
print(df_long)

from scipy.stats import mannwhitneyu
import statannotations.Annotator

# Set style
sns.set(style="whitegrid", font_scale=1.1)

# Create grouped violin plot
plt.figure(figsize=(6, 5))
ax = sns.violinplot(
    data=df_long,
    x="comparison",
    y="log10_length",
    hue="species",
    dodge=True,
    inner="quartile",
    palette="Set2",
    width=0.5
)

# Define pairs for annotation (must match violin grouping)
pairs = [
    (("OE_vs_CV", "C. virginica"), ("OE_vs_CV", "Other oysters")),
    (("MG_vs_CV", "C. virginica"), ("MG_vs_CV", "Other oysters")),
    (("MA_vs_CV", "C. virginica"), ("MA_vs_CV", "Other oysters"))
]

# Initialize the annotator
annotator = Annotator(ax, pairs, data=df_long, x="comparison", y="log10_length", hue="species", dodge=True)

# Add statistical test (Mann–Whitney U) and annotate with stars
annotator.configure(test="Mann-Whitney", text_format="star", loc="outside", comparisons_correction=None)
annotator.apply_and_annotate()

# Finalize plot
plt.ylabel("log10(Synteny Block Length)")
plt.xlabel("Comparison Genome")
ax.grid(False)
sns.despine(top=True,right=True)
ax.set_ylim(4,8.5)
ax.set_xticklabels(["O. edulis","M. gigas","M. angulata"])
ax.legend(title="Genome",loc="upper right",frameon=False)
plt.tight_layout()
plt.savefig("violin_synteny_block_lengths.png",dpi=600)



# Define chromosome lengths
chromosome_lengths = {
    "cf1": 64991568,
    "cf2": 59045500,
    "cf3": 58272962,
    "cf4": 56786176,
    "cf5": 55461147,
    "cf6": 50076034,
    "cf7": 49270460,
    "cf8": 48947930,
    "cf9": 48412105,
    "cf10": 37762054
}

# Load your synteny block table
df = pd.read_csv("synteny_block_lengths.tsv", sep="\t")  # Update with your actual file


chrom_order = [f"cf{i}" for i in range(1, 11)]  # ["cf1", ..., "cf10"]
chrom_labels = [f"chr{i}" for i in range(1, 11)]  # ["chr1", ..., "chr10"]

comparison_label_map = {
    "OE_vs_CV": "O. edulis",
    "MG_vs_CV": "M. gigas",
    "MA_vs_CV": "M. angulata"
}

# Group by comparison and chromosome
synteny_sum = (
    df.groupby(["comparison", "cv_chromosome"])["block_length_cv"]
    .sum()
    .reset_index()
    .rename(columns={"block_length_cv": "total_syntenic_length"})
)

# Add chromosome lengths
synteny_sum["chrom_length"] = synteny_sum["cv_chromosome"].map(chromosome_lengths)

# Calculate percentage
synteny_sum["synteny_percent"] = 100 * synteny_sum["total_syntenic_length"] / synteny_sum["chrom_length"]
# Create a new column with cleaned legend labels
synteny_sum["comparison_clean"] = synteny_sum["comparison"].map(comparison_label_map)

# Make sure chromosome order is categorical for sorting
synteny_sum["cv_chromosome"] = pd.Categorical(
    synteny_sum["cv_chromosome"],
    categories=chrom_order,
    ordered=True
)

# Re-sort dataframe
synteny_sum = synteny_sum.sort_values("cv_chromosome")
hue_order = ["O. edulis","M. gigas", "M. angulata"]

print(synteny_sum)

# Plot
plt.figure(figsize=(6, 5))
ax = sns.barplot(
    data=synteny_sum,
    x="cv_chromosome",
    y="synteny_percent",
    hue="comparison_clean",
    palette="Set2",
    hue_order=hue_order
)

# Rename x-ticks from cf1 → chr1, cf2 → chr2, ...
ax.set_xticklabels(chrom_labels)

plt.ylabel("Synteny Coverage (%)")
plt.xlabel("C. virginica Chromosome")
plt.legend(title="Genome",frameon=False,loc='upper right')
ax.set_ylim(0,140)
yticks = ax.get_yticks()
ax.set_yticks([tick for tick in yticks if tick not in [120, 140]])
ax.grid(False)
sns.despine(top=True,right=True)
plt.tight_layout()
plt.savefig("Synteny_percentage_by_chromosome.png",dpi=300)

# Pivot to wide format
pivot_df = synteny_sum.pivot(index="cv_chromosome", columns="comparison_clean", values="synteny_percent")

from scipy.stats import wilcoxon

stat, pval = wilcoxon(pivot_df["O. edulis"], pivot_df[["M. angulata", "M. gigas"]].mean(axis=1), alternative="less")
print(f"Wilcoxon signed-rank: W = {stat:.2f}, p = {pval:.4e}")

# Label chromosome group
synteny_sum["chr_group"] = synteny_sum["cv_chromosome"].apply(
    lambda x: "chr3/10" if x in ["cf3", "cf10"] else "other"
)

# Extract values
group_a = synteny_sum[synteny_sum["chr_group"] == "chr3/10"]["synteny_percent"]
group_b = synteny_sum[synteny_sum["chr_group"] == "other"]["synteny_percent"]

stat2, pval2 = mannwhitneyu(group_a, group_b, alternative="less")

print("Test: chr3/10 vs other chromosomes")
print(f"U = {stat2:.2f}, p = {pval2:.4e}")

# # === Grouped summary statistics ===
# summary_stats = df_long.groupby(["comparison", "species"]).agg(
#     count=("block_length", "count"),
#     mean=("block_length", "mean"),
#     median=("block_length", "median"),
#     std=("block_length", "std"),
#     min=("block_length", "min"),
#     max=("block_length", "max")
# ).reset_index()

# # Display in console
# print("\n=== Syntenic Block Length Summary Stats ===")
# print(summary_stats)

# # Save to TSV
# summary_stats.to_csv("synteny_block_stats_by_group.tsv", sep="\t", index=False)

