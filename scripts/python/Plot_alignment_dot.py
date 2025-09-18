import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")


# === Settings ===
coords_file = "cv_vs_cv_new_all.filtered.coords"  # your coords file
MIN_ALIGN_LEN_MB = 0.01   # Only keep alignments >50 kb
MIN_IDENTITY = 90         # Only keep alignments >95% identity
# === Load and preprocess ===

columns = [
    "Ref_start", "Ref_end", "Qry_start", "Qry_end",
    "Ref_len", "Qry_len", "Identity", "Ref_chr_num",
    "Strand", "Ref_chr", "Qry_chr"
]
df = pd.read_csv(coords_file, delim_whitespace=True, header=None, names=columns)

# Compute alignment length (Mb)
df["Align_len"] = abs(df["Ref_end"] - df["Ref_start"]) / 1e6  # in Mb

df = df[(df["Align_len"] > MIN_ALIGN_LEN_MB) & (df["Identity"] > MIN_IDENTITY)]

# Create pivot table: rows = old reference, columns = new assembly
pivot = df.groupby(["Ref_chr", "Qry_chr"])["Align_len"].sum().unstack().fillna(0)

# Sort axes
def chr_key(x): return int(x.replace("chr", "")) if x.startswith("chr") else 999
def ncbi_key(x):
    try:
        return int(x.split("_")[1].split(".")[0])
    except (IndexError, ValueError):
        return 999

pivot = pivot.loc[sorted(pivot.index, key=ncbi_key), sorted(pivot.columns, key=chr_key)]
pivot = pivot.iloc[::-1]  # Flip y-axis (increasing chr1 → chr10 from bottom up)

# === Plotting ===
sns.set(style="whitegrid", font_scale=1.2)
fig, ax = plt.subplots(figsize=(12, 12))

heatmap = sns.heatmap(
    pivot,
    cmap="YlGnBu",
    linewidths=0.5,
    linecolor='gray',
    annot=True,
    fmt=".1f",
    annot_kws={"fontsize": 10},
    cbar_kws={"label": "Alignment Length (Mb)","shrink":0.3,"pad":0.01},
    square=True,
    ax=ax
)

# Axis labels
ax.set_xlabel("New Assembly", fontsize=14)
ax.set_ylabel("C.virginica_3.0", fontsize=14)

# Tick formatting
plt.xticks(rotation=0, ha="right", fontsize=11)
plt.yticks(rotation=0, fontsize=11)

# Save high-res
plt.tight_layout()
plt.savefig("chr_alignment_heatmap.png", dpi=600)


# === Step 1: Total alignment lengths per pair ===
pairwise = df.groupby(["Ref_chr", "Qry_chr"])["Align_len"].sum().reset_index()

# === Step 2: Total alignment per Ref_chr ===
ref_total = df.groupby("Ref_chr")["Align_len"].sum().rename("Ref_total_len").reset_index()

# === Step 3: Find top query hit per Ref_chr ===
top_qry_per_ref = pairwise.sort_values("Align_len", ascending=False).groupby("Ref_chr").first().reset_index()
# print(top_qry_per_ref)

# Merge with total to compute fraction
ref_merged = top_qry_per_ref.merge(ref_total, on="Ref_chr")
ref_merged["TopHit_fraction"] = ref_merged["Align_len"] / ref_merged["Ref_total_len"]

# === Step 4: Repeat for Query chromosomes ===
qry_total = df.groupby("Qry_chr")["Align_len"].sum().rename("Qry_total_len").reset_index()
top_ref_per_qry = pairwise.sort_values("Align_len", ascending=False).groupby("Qry_chr").first().reset_index()
qry_merged = top_ref_per_qry.merge(qry_total, on="Qry_chr")
qry_merged["TopHit_fraction"] = qry_merged["Align_len"] / qry_merged["Qry_total_len"]

# print(qry_merged["TopHit_fraction"])
# print(ref_merged["TopHit_fraction"])


###Dot plot
# Minimum alignment filter
ref_fai_file = "Cv3.0.fa.fai"
qry_fai_file = "primary_dedup_chr_masked_hp_sealed.fa.fai"


# === 1. Load alignment data ===
columns = [
    "Ref_start", "Ref_end", "Qry_start", "Qry_end",
    "Ref_len", "Qry_len", "Identity", "Ref_chr_num",
    "Strand", "Ref_chr", "Qry_chr"
]
df = pd.read_csv(coords_file, delim_whitespace=True, header=None, names=columns)

# Compute alignment length in Mb
df["Align_len"] = abs(df["Ref_end"] - df["Ref_start"]) / 1e6

# Filter for quality
df = df[(df["Align_len"] > MIN_ALIGN_LEN_MB) & (df["Identity"] > MIN_IDENTITY)]
print(f"Filtered alignments retained: {len(df)}")

# === 2. Load .fai files for chromosome sizes ===
ref_fai = pd.read_csv(ref_fai_file, sep="\t", header=None,usecols=[0,1], names=["chrom","length"])
qry_fai = pd.read_csv(qry_fai_file, sep="\t", header=None,usecols=[0,1], names=["chrom","length"])
# Manual renaming map
qry_rename_map = {
    "HiC_scaffold_1": "chr1",
    "HiC_scaffold_2": "chr2",
    "HiC_scaffold_3": "chr3",
    "HiC_scaffold_4": "chr4",
    "HiC_scaffold_5": "chr5",
    "HiC_scaffold_6": "chr6",
    "HiC_scaffold_7": "chr7",
    "HiC_scaffold_8": "chr8",
    "HiC_scaffold_9": "chr9",
    "HiC_scaffold_10": "chr10"
}

# Apply it
qry_fai["chrom"] = qry_fai["chrom"].map(qry_rename_map).fillna(qry_fai["chrom"])

ref_fai["length_mb"] = ref_fai["length"] / 1e6
qry_fai["length_mb"] = qry_fai["length"] / 1e6

# Filter to only chromosomes actually in alignments
ref_fai = ref_fai[ref_fai["chrom"].isin(df["Ref_chr"].unique())].copy()
qry_fai = qry_fai[qry_fai["chrom"].isin(df["Qry_chr"].unique())].copy()

# Sort chromosomes
def parse_ref_chr(x):
    try:
        return int(x.split("_")[1].split(".")[0])
    except:
        return 999

def parse_qry_chr(x):
    try:
        return int(x.replace("chr",""))
    except:
        return 999

ref_fai = ref_fai.sort_values(by="chrom", key=lambda s: s.map(parse_ref_chr))
qry_fai = qry_fai.sort_values(by="chrom", key=lambda s: s.map(parse_qry_chr))

# Compute cumulative offsets
ref_fai["offset"] = ref_fai["length_mb"].cumsum().shift(fill_value=0)
qry_fai["offset"] = qry_fai["length_mb"].cumsum().shift(fill_value=0)

# Map offsets onto alignment dataframe
ref_offset_map = ref_fai.set_index("chrom")["offset"].to_dict()
qry_offset_map = qry_fai.set_index("chrom")["offset"].to_dict()

df["Ref_offset"] = df["Ref_chr"].map(ref_offset_map)
df["Qry_offset"] = df["Qry_chr"].map(qry_offset_map)

# Linearized genome-wide coordinates
df["Ref_linear"] = (df["Ref_start"] / 1e6) + df["Ref_offset"]
df["Qry_linear"] = (df["Qry_start"] / 1e6) + df["Qry_offset"]

# === 3. Plotting ===
plt.figure(figsize=(12,12))

ax = plt.gca()

# Colors for strand orientation
strand_colors = {1: "blue", -1: "red"}
for strand in [-1, 1]:
    subset = df[df["Strand"] == strand]
    ax.scatter(
        subset["Qry_linear"],
        subset["Ref_linear"],
        s=0.3,
        c=strand_colors[strand],
        label=f"Strand {strand}",
        alpha=0.4
    )

# Chromosome grid lines
for x in qry_fai["offset"]:
    ax.axvline(x, color='black', lw=1, linestyle='--', alpha=0.7)
for y in ref_fai["offset"]:
    ax.axhline(y, color='black', lw=1, linestyle='--', alpha=0.7)


# Axis tick labels at midpoints
ref_midpoints = ref_fai["offset"] + ref_fai["length_mb"] / 2
qry_midpoints = qry_fai["offset"] + qry_fai["length_mb"] / 2


ax.set_xticks(qry_midpoints)
ax.set_xticklabels(qry_fai["chrom"],fontsize=13)
ax.set_yticks(ref_midpoints)
ax.set_yticklabels(ref_fai["chrom"],fontsize=13)


# Axis limits to span full genome
ax.set_xlim(0, qry_fai["offset"].iloc[-1] + qry_fai["length_mb"].iloc[-1])
ax.set_ylim(0, ref_fai["offset"].iloc[-1] + ref_fai["length_mb"].iloc[-1])
ax.margins(0)


# Labels, legend, title
ax.set_xlabel("New Assembly", fontsize=13)
ax.set_ylabel("C.virginica_3.0", fontsize=13)
ax.legend(title="Alignment Strand", labels=["Reverse","Forward"], fontsize=12, title_fontsize=14, markerscale=6)
ax.grid(False)


# Save outputs
plt.tight_layout(pad=0)
plt.savefig("genomewide_dotplot_matrix_publication.png", dpi=600, bbox_inches='tight')


# === 4. Alignment summary statistics ===
total_alignment_mb = df["Align_len"].sum()
print(f"Total alignment length: {total_alignment_mb:.2f} Mb")


ref_total_size_mb = ref_fai["length_mb"].sum()
ref_coverage_fraction = total_alignment_mb / ref_total_size_mb
print(f"Fraction of reference covered: {ref_coverage_fraction:.2%}")


qry_total_size_mb = qry_fai["length_mb"].sum()
qry_coverage_fraction = total_alignment_mb / qry_total_size_mb
print(f"Fraction of query covered: {qry_coverage_fraction:.2%}")


print("\nAlignment identity stats:")
print(df["Identity"].describe())


print("\nStrand orientation fractions:")
strand_counts = df["Strand"].value_counts(normalize=True)
print(strand_counts)

