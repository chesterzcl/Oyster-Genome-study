import pandas as pd
import os
from statannotations.Annotator import Annotator


os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final")

# ================================
# CONFIG
# ================================
ORTHOGROUP_FILE = "Orthogroups.tsv"
BREAKEND_OVERLAP_FILE = "breakend_flank_gene_overlaps.csv"
GENE_FUNCTION_FILE = "genes_function.bed"

OUTPUT_GENE_ANNOT = "all_genes_orthogroup_categories.csv"
OUTPUT_BREAKEND_ANNOT = "breakend_flank_genes_with_orthogroup_category.csv"

# ================================
# 1. Load Orthogroups.tsv
# ================================
df_ortho = pd.read_csv(ORTHOGROUP_FILE, sep="\t")

# Columns: Orthogroup, CV, MG
cv_col = [col for col in df_ortho.columns if 'CV' in col or 'cv' in col][0]
mg_col = [col for col in df_ortho.columns if 'MG' in col or 'mg' in col][0]

# Parse gene lists
def parse_genes(cell):
    if pd.isna(cell): return []
    return [x.strip() for x in cell.split(",") if x.strip()]

df_ortho["CV_genes"] = df_ortho[cv_col].apply(parse_genes)
df_ortho["MG_genes"] = df_ortho[mg_col].apply(parse_genes)

# Count copies
df_ortho["CV_count"] = df_ortho["CV_genes"].apply(len)
df_ortho["MG_count"] = df_ortho["MG_genes"].apply(len)

# ================================
# 2. Assign orthogroup categories
# ================================
def classify_row(row):
    cv = row["CV_count"]
    mg = row["MG_count"]
    
    if cv == 0 and mg == 0:
        return "Empty"
    elif cv > 0 and mg == 0:
        return "CV_specific"
    elif cv == 0 and mg > 0:
        return "MG_specific"
    elif cv == 1 and mg == 1:
        return "1:1"
    elif cv > 1 and mg == 1:
        return "CV_expansion"
    elif cv == 1 and mg > 1:
        return "MG_expansion"
    elif cv > 1 and mg > 1:
        return "Many:Many"
    else:
        return "Other"

df_ortho["category"] = df_ortho.apply(classify_row, axis=1)
print(df_ortho["category"].value_counts())

# ================================
# 3. Build gene -> category mapping (from orthogroups)
# ================================
records = []
for _, row in df_ortho.iterrows():
    if isinstance(row["CV_genes"], list):
        for gene in row["CV_genes"]:
            gene = gene.split('|')[1] if '|' in gene else gene
            if isinstance(gene, str) and gene.strip():
                records.append([gene, row["Orthogroup"], row["category"], "CV"])
    if isinstance(row["MG_genes"], list):
        for gene in row["MG_genes"]:
            gene = gene.split('|')[1] if '|' in gene else gene
            if isinstance(gene, str) and gene.strip():
                records.append([gene, row["Orthogroup"], row["category"], "MG"])

gene_ortho_df = pd.DataFrame(records, columns=["gene_id", "orthogroup", "category", "species"])

# ================================
# 4. Add "Unassigned" category for CV genes without orthogroup
# ================================
cv_gene_annot = pd.read_csv(
    GENE_FUNCTION_FILE, sep="\t", header=None,
    names=["chrom", "start", "end", "gene_id", "gene_name", "strand"]
)

# Clean gene_id
cv_gene_annot["gene_id"] = cv_gene_annot["gene_id"].astype(str).str.strip()

all_cv_genes = set(cv_gene_annot["gene_id"])
mapped_cv_genes = set(gene_ortho_df[gene_ortho_df["species"]=="CV"]["gene_id"])

unassigned_genes = all_cv_genes - mapped_cv_genes
print(f"Number of unassigned CV genes: {len(unassigned_genes)}")

unassigned_records = [
    [gene, None, "Unassigned", "CV"] for gene in unassigned_genes
]
unassigned_df = pd.DataFrame(
    unassigned_records,
    columns=["gene_id", "orthogroup", "category", "species"]
)

gene_ortho_df = pd.concat([gene_ortho_df, unassigned_df], ignore_index=True)
gene_ortho_df.to_csv(OUTPUT_GENE_ANNOT, index=False)
print(f"*** Saved {OUTPUT_GENE_ANNOT} ***")

# ================================
# 5. Merge with breakend flank genes
# ================================
flank_df = pd.read_csv(BREAKEND_OVERLAP_FILE)
flank_df["gene_id"] = flank_df["gene_id"].astype(str).str.strip()

merged_df = flank_df.merge(
    gene_ortho_df[gene_ortho_df["species"]=="CV"],
    on="gene_id",
    how="left"
)

# Fill any missing category
merged_df["category"] = merged_df["category"].fillna("Unassigned")
merged_df.to_csv(OUTPUT_BREAKEND_ANNOT, index=False)
print(f"*** Saved {OUTPUT_BREAKEND_ANNOT} ***")


import matplotlib.pyplot as plt
import seaborn as sns

# ================================
# 6. Count categories at breakends
# ================================
breakend_annot = pd.read_csv(OUTPUT_BREAKEND_ANNOT)

# Count the number of genes in each category
category_counts = (
    breakend_annot
    .groupby("category")
    .size()
    .reset_index(name="count")
    .sort_values("count", ascending=False)
)

# Save
category_counts.to_csv("breakend_category_counts.csv", index=False)
print("*** Saved breakend_category_counts.csv ***")

# ================================
# 7. Count genes per flank, per category
# ================================
per_flank_counts = (
    breakend_annot
    .groupby(["flank_name", "event_type", "category"])
    .size()
    .reset_index(name="gene_count")
)

per_flank_counts.to_csv("breakend_per_flank_gene_counts.csv", index=False)
print("*** Saved breakend_per_flank_gene_counts.csv ***")


# ================================
# 8. Plot grouped boxplot
# ================================
import seaborn as sns
import matplotlib.pyplot as plt

# Remove MG_expansion category
plot_df = per_flank_counts[
    (per_flank_counts["category"] != "MG_expansion") &
    (per_flank_counts["category"] != "Unassigned")
]

plt.figure(figsize=(10,6))
ax = sns.boxplot(
    data=plot_df,
    x="category",
    y="gene_count",
    hue="event_type",
    showfliers=False,
    width=0.6,
    hue_order=["inversion","insertion_cv","insertion_mg","insertion_mg_foreign", "complex_sv" ]
    # palette="Set2"
)
# Custom labels
ax.set_xticklabels([
    "Single-copy",
    "CV expansion",
    "CV-specific",
    "Reciprocal expansion"
])
plt.xlabel("Orthogroup Category")
plt.ylabel("Number of Genes per 100kb Flank Window")


# Customize legend labels
handles, labels = ax.get_legend_handles_labels()

# Example: Replace labels with custom ones
custom_labels = ['Inversion','Deletion','Local Insertion','Foreign Insertion', 'Complex SV']
ax.legend(handles=handles, labels=custom_labels, title=None, loc='upper right')
# ax.legend(handles=handles, title="SV Type", loc='upper right')

plt.tight_layout()
plt.savefig("breakend_category_boxplot_by_event_type.png")
plt.close()

print("*** Saved breakend_category_boxplot_by_event_type.png ***")



INPUT_GAP_CSV="gap_records_analysis.csv"
# === Step 1: Load input files ===
gap_df = pd.read_csv(INPUT_GAP_CSV)

# Only keep insertions
# insertion_df = gap_df[gap_df["event_type"].isin(["insertion_cv", "insertion_mg"])].copy()
insertion_df = gap_df[(gap_df["event_type"] != "trivial") ].copy()
# Replace cf_chr with HiC scaffold names
cf_to_hic = {f"cf{i}": f"HiC_scaffold_{i}" for i in range(1, 11)}    
insertion_df["cf_chr"] = insertion_df["cf_chr"].map(cf_to_hic)
cv_gene_annot["gene_id"] = cv_gene_annot["gene_id"].astype(str).str.strip()

# Merge with orthogroup category
cv_ortho = gene_ortho_df[gene_ortho_df["species"] == "CV"]
gene_annot = cv_gene_annot.merge(cv_ortho, on="gene_id", how="left")
gene_annot["category"] = gene_annot["category"].fillna("Unassigned")
print(list(cv_gene_annot))
print(len(gene_annot))
# === Step 2: Scan insertion gaps for overlapping genes ===
records = []

for idx, row in insertion_df.iterrows():
    chrom = row["cf_chr"]
    start = min(row["cf_pos1"], row["cf_pos2"])
    end = max(row["cf_pos1"], row["cf_pos2"])
    event_type = row["event_type"]

    gap_id = f"{chrom}:{start}-{end}"
    gap_len=abs(end-start)

    overlapping = gene_annot[
        (gene_annot["chrom"] == chrom) &
        (gene_annot["end"] > start) &
        (gene_annot["start"] < end)
    ]

    for cat, count in overlapping["category"].value_counts().items():
        records.append([gap_id, event_type, cat, count,gap_len])

# === Step 3: Save result ===
gap_cat_df = pd.DataFrame(records, columns=["gap_id", "event_type", "category", "gene_count","gap_length"])
gap_cat_df.to_csv("gap_insertion_orthogroup_distribution.csv", index=False)
print("*** Saved gap_insertion_orthogroup_distribution.csv ***")

# Optional: keep major categories only
plot_df = gap_cat_df[gap_cat_df["category"].isin(["1:1","CV_specific", "CV_expansion", "Many:Many"])].copy()
# print(plot_df)
# Rename for clean legend and x-axis
event_rename = {
    "insertion_cv": "Deletion",
    "insertion_mg": "Local Insertion",
    # "insertion_mg_foreign":"Foreign Insertion",
    # "complex_sv":"Complex",
    "inversion":"Inversion"
}
plot_df["event_type"] = plot_df["event_type"].map(event_rename)
# print(plot_df)
category_order = ["1:1","CV_specific", "CV_expansion", "Many:Many"]
plot_df["category"] = pd.Categorical(plot_df["category"], categories=category_order, ordered=True)
plot_df["gene_density"]=plot_df["gene_count"]/plot_df["gap_length"]
# Plot
plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=plot_df,
    x="category",
    # y="gene_density",
    y="gene_count",
    hue="event_type",
    showfliers=False,
    width=0.6
)

# # Define the pairs to compare for each category
# categories = plot_df['category'].unique()
# event_types = plot_df['event_type'].unique()

# # Generate all pairwise comparisons within each category
# pairs = []
# for cat in categories:
#     events_in_cat = plot_df[plot_df['category'] == cat]['event_type'].unique()
#     if len(events_in_cat) >= 2:
#         for i in range(len(events_in_cat)):
#             for j in range(i+1, len(events_in_cat)):
#                 pairs.append(((cat, events_in_cat[i]), (cat, events_in_cat[j])))

# # Add statistical annotations
# annotator = Annotator(ax, pairs, data=plot_df, x="category", y="gene_count", hue="event_type")
# annotator.configure(test="Mann-Whitney", text_format="star", loc="inside", comparisons_correction="bonferroni")
# annotator.apply_and_annotate()


# Styling
plt.ylabel("Number of Genes per Insertion")
ax.set_xticklabels([
    "Single-copy",
    "CV expansion",
    "CV-specific",
    "Reciprocal expansion"
])
plt.xlabel("Orthogroup Category")
# plt.title("Distribution of Orthogroup Types per Insertion Event")
plt.legend(title="Insertion Type", loc='upper right')
handles, labels = ax.get_legend_handles_labels()
custom_labels = ["CV Insertion","MG Insertion"]
# ax.legend(handles=handles, labels=custom_labels, title=None, loc='upper right')
ax.legend(handles=handles, title=None, loc='upper right')

plt.tight_layout()
plt.savefig("orthogroup_category_boxplot_per_insertion.png", dpi=300)
plt.close()

print("*** Saved orthogroup_category_boxplot_per_insertion.png ***")

####=================================================================================
# Step 1: Define mapping function
def map_to_expansion_only(cat):
    if cat in ["CV_specific", "CV_expansion", "MG_expansion", "Many:Many"]:
        return "Expansion/Specific"
    elif cat == "1:1":
        return "1:1"
    else:
        return None  # Will be excluded

# Step 2: Apply the mapping
filtered_df = gap_cat_df.copy()
filtered_df["expansion_status"] = filtered_df["category"].apply(map_to_expansion_only)

# Step 3: Drop rows that are not classified (i.e., not expansion or 1:1)
filtered_df = filtered_df.dropna(subset=["expansion_status"])

# Step 4: Aggregate gene counts by new category
summary_df = (
    filtered_df
    .groupby(["gap_id", "event_type", "expansion_status", "gap_length"], as_index=False)
    .agg({"gene_count": "sum"})
)

# (Optional) Save the result
summary_df.to_csv("gap_insertion_expansion_vs_1to1.csv", index=False)
print("*** Saved gap_insertion_expansion_vs_1to1.csv ***")


# Pivot to wide format: one row per gap, columns for expansion and 1:1 gene counts
pivot_df = (
    summary_df
    .pivot_table(
        index=["gap_id", "event_type", "gap_length"],
        columns="expansion_status",
        values="gene_count",
        aggfunc="sum",
        fill_value=0  # ← ensures missing combinations are filled with 0
    )
    .reset_index()
)

# Clean column names
pivot_df.columns.name = None
pivot_df = pivot_df.rename(columns={
    "Expansion/Specific": "expansion_genes",
    "1:1": "one2one_genes"
})

# Rename columns for clarity
pivot_df.columns.name = None  # remove pandas category name
pivot_df = pivot_df.rename(columns={
    "Expansion/Specific": "expansion_genes",
    "1:1": "one2one_genes"
})

from scipy.stats import mannwhitneyu

sv_types = pivot_df["event_type"].unique()

for sv in sv_types:
    sub = pivot_df[pivot_df["event_type"] == sv]
    
    # Extract expansion and 1:1 gene counts
    expansion_values = sub["expansion_genes"]
    one2one_values = sub["one2one_genes"]
    
    # Filter out rows where both counts are 0
    valid = sub[(sub["expansion_genes"] > 0) | (sub["one2one_genes"] > 0)]

    if len(valid) < 3:
        print(f"{sv}: Not enough data for test.")
        continue

    # Re-extract filtered values
    expansion_values = valid["expansion_genes"]
    one2one_values = valid["one2one_genes"]

    try:
        stat, p = mannwhitneyu(expansion_values, one2one_values, alternative='two-sided')
        print(f"{sv}: Mann–Whitney U test p = {p:.4g} (n_exp = {len(expansion_values)}, n_1to1 = {len(one2one_values)})")
    except ValueError as e:
        print(f"{sv}: Test failed — {str(e)}")



# Melt the wide pivot_df into long format
melted_df = pivot_df.melt(
    id_vars=["gap_id", "event_type"],
    value_vars=["expansion_genes", "one2one_genes"],
    var_name="Orthogroup Type",
    value_name="Gene Count"
)

# Make Orthogroup Type values nicer
melted_df["Orthogroup Type"] = melted_df["Orthogroup Type"].replace({
    "expansion_genes": "Expansion/Specific",
    "one2one_genes": "1:1"
})

sv_keep = ["insertion_cv", "insertion_mg", "insertion_mg_foreign","inversion","complex_sv"]
plot_df = melted_df[melted_df["event_type"].isin(sv_keep)].copy()

print("X values:", plot_df["event_type"].unique())
print("Hue values:", plot_df["Orthogroup Type"].unique())

x_order = ["complex_sv", "insertion_cv", "insertion_mg", "insertion_mg_foreign", "inversion"]
hue_order = ["Expansion/Specific", "1:1"]

plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=plot_df,
    x="event_type",
    y="Gene Count",
    hue="Orthogroup Type",
    order=x_order,
    hue_order=hue_order,
    showfliers=False,
    width=0.6
)

# Styling
ax.set_xlabel("SV Gap Type", fontsize=14)
ax.set_ylabel("Gene Count per Gap", fontsize=14)
# ax.set_title("Gene Expansion vs 1:1 per SV Type")
# Rename x-tick labels (example)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles,title="Orthogroup Category", labels=["Expanded/Lineage-specific", "One-to-one"])

from statannotations.Annotator import Annotator

pairs = [
    (("complex_sv", "Expansion/Specific"), ("complex_sv", "1:1")),
    (("insertion_cv", "Expansion/Specific"), ("insertion_cv", "1:1")),
    (("insertion_mg", "Expansion/Specific"), ("insertion_mg", "1:1")),
    (("insertion_mg_foreign", "Expansion/Specific"), ("insertion_mg_foreign", "1:1")),
    (("inversion", "Expansion/Specific"), ("inversion", "1:1")),
]

# Corresponding p-values from your test results
manual_pvals = [
    9.928e-05,
    2.115e-27,
    0.06883,
    0.008862,
    1.678e-05
]

# Apply manual annotation
annotator = Annotator(ax, pairs, data=plot_df, x="event_type", y="Gene Count", hue="Orthogroup Type",order=x_order, hue_order=hue_order)
annotator.configure(test=None, text_format='star')  # test=None disables internal testing
annotator.set_pvalues_and_annotate(pvalues=manual_pvals)

ax.set_xticklabels([
    "Complex SV",
    "Deletion",
    "Local Insertion",
    "Foreign Insertion",
    "Inversion"
])
plt.tight_layout()
plt.savefig("expansion_vs_1to1_by_svtype_boxplot.png", dpi=600)
plt.close()
####=================================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# === Input ===
fai_file = "primary_dedup_chr_masked_hp_sealed.fa.fai"
chrom_sizes = pd.read_csv(fai_file, sep="\t", usecols=[0, 1], names=["chrom", "length"], header=None)
fs=20
# === Color Maps ===
category_colors = {
    "1:1": "#4daf4a",
    "CV_specific": "#e41a1c",
    "CV_expansion": "#377eb8",
    "Many:Many": "#984ea3",
    "Unassigned": "gray"
}
sv_colors = {
    "inversion": "#fdb863",           # orange-ish yellow
    "insertion_cv": "#b2abd2",        # soft lavender
    "insertion_mg": "#5e3c99",        # dark purple
    "insertion_mg_foreign": "#e66101",# orange
    "complex_sv": "#1b7837"           # dark green
}

window_size = 500_000

# === Bin genes by orthogroup category ===
gene_annot = gene_annot[gene_annot["end"] - gene_annot["start"] > 400]
bins = []
for chrom, length in chrom_sizes.values:
    for start in range(0, length, window_size):
        end = start + window_size
        mid = (start + end) // 2

        genes_in_bin = gene_annot[
            (gene_annot["chrom"] == chrom) &
            (gene_annot["start"] < end) &
            (gene_annot["end"] > start)
        ]

        record = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "mid": mid
        }
        for cat in category_colors:
            record[cat] = (genes_in_bin["category"] == cat).sum()

        bins.append(record)

bin_df = pd.DataFrame(bins)
# Sum all category counts in each bin
bin_df["total_count"] = bin_df[list(category_colors.keys())].sum(axis=1)
ymax = bin_df["total_count"].max()
# =============== Summarizing heterozygosity ========================
# Add per-bin ratio
bin_df['1:1_ratio'] = bin_df['1:1'] / bin_df['total_count']

# Per-chromosome summary
chrom_summary = (
    bin_df.groupby('chrom')['1:1_ratio']
    .agg(['mean', 'std', 'min', 'max', 'median', 'count'])
    .reset_index()
    .rename(columns={
        'mean': 'mean_1:1_ratio',
        'std': 'std_1:1_ratio',
        'min': 'min_1:1_ratio',
        'max': 'max_1:1_ratio',
        'median': 'median_1:1_ratio',
        'count': 'num_bins'
    })
)

# Genome-wide summary
genome_summary = bin_df['1:1_ratio'].agg(['mean', 'std', 'min', 'max', 'median', 'count']).to_frame().T
genome_summary['chrom'] = 'Genome'
genome_summary = genome_summary[[
    'chrom', 'mean', 'std', 'min', 'max', 'median', 'count'
]].rename(columns={
    'mean': 'mean_1:1_ratio',
    'std': 'std_1:1_ratio',
    'min': 'min_1:1_ratio',
    'max': 'max_1:1_ratio',
    'median': 'median_1:1_ratio',
    'count': 'num_bins'
})

# Combine and calculate Coefficient of Variation
summary_with_genome = pd.concat([chrom_summary, genome_summary], ignore_index=True)
summary_with_genome['cv_1:1_ratio'] = summary_with_genome['std_1:1_ratio'] / summary_with_genome['mean_1:1_ratio']

summary_with_genome.to_csv("1to1_ratio_summary_per_chromosome.csv", index=False)




# ===================================================================
# === Plotting ===

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

# Determine max chromosome length for shared X axis
max_length = chrom_sizes["length"].max()

# Map cf1 → HiC_scaffold_1, etc.
chrom_map = {f"cf{i}": f"HiC_scaffold_{i}" for i in range(1, 11)}
# Build mapping from "HiC_scaffold_1" → "chr1"
chrom_label_map = {
    f"HiC_scaffold_{i}": f"chr{i}" for i in range(1, 11)
}
# Apply mapping
gap_df["cf_chr"] = gap_df["cf_chr"].map(chrom_map)

# === Plotting ===
fig = plt.figure(figsize=(16, 3 * len(chrom_sizes)))
gs = gridspec.GridSpec(len(chrom_sizes)*2, 1, height_ratios=[1.2,0.5] * len(chrom_sizes),hspace=0.1)

for i, (chrom, length) in enumerate(chrom_sizes.values):
    ax_bar = fig.add_subplot(gs[i*2])
    ax_chr = fig.add_subplot(gs[i*2 + 1])

    # === Top track: orthogroup stacked bar ===
    sub_bin = bin_df[bin_df["chrom"] == chrom]
    base_y = 0
    for cat in category_colors:
        height = sub_bin[cat]
        ax_bar.bar(sub_bin["mid"], height, width=window_size * 0.9,
                   bottom=base_y, color=category_colors[cat],
                   label=cat if i == 0 else None, alpha=0.8, linewidth=0)
        base_y += height

    ax_bar.set_xlim(0, max_length)
    ax_bar.set_xticks([])
    ax_bar.set_yticks([])
    ax_bar.set_ylim(0, ymax * 1.05)
    # ax_bar.set_title(chrom, loc='left', fontsize=12)
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.spines['left'].set_visible(False)
    ax_bar.spines['bottom'].set_visible(False)

    # === Bottom track: chromosome bar + SV gaps ===
    ax_chr.add_patch(patches.Rectangle((0, 0.1), length, 0.2, color='gray'))

    sub_gaps = gap_df[gap_df["cf_chr"] == chrom]
    for _, row in sub_gaps.iterrows():
        color = sv_colors.get(row["event_type"], "gray")
        ax_chr.add_patch(patches.Rectangle(
            (row["cf_pos1"], 0.1), row["cf_pos2"] - row["cf_pos1"], 0.2,
            color=color, alpha=1
        ))

    ax_chr.set_xlim(0, max_length)
    ax_chr.set_ylim(0, 0.5)
    ax_chr.set_yticks([])
    if i == len(chrom_sizes) - 1:
        ax_chr.set_xlabel("Genomic position (Mb)",fontsize=fs)
        ax_chr.xaxis.set_major_locator(ticker.MultipleLocator(10_000_000))
        ax_chr.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.0f}"))
    else:
        ax_chr.set_xticks([])

    label = chrom_label_map.get(chrom, chrom)
    ax_chr.text(-0.01 * max_length, 0.25, label, ha='right', va='center', fontsize=fs, transform=ax_chr.transData)
    
    ax_chr.spines['top'].set_visible(False)
    ax_chr.spines['right'].set_visible(False)
    ax_chr.spines['left'].set_visible(False)
    if i != len(chrom_sizes) - 1:
        ax_chr.spines['bottom'].set_visible(False)
        ax_chr.tick_params(labelsize=18)

# === Legend ===
import matplotlib.patches as mpatches
from matplotlib.patches import Patch

orthogroup_labels = {
    "1:1": "One-to-One",
    "CV_specific": "CV-Specific",
    "CV_expansion": "CV Expansion",
    "Many:Many": "Many-to-Many",
    "Unassigned": "Unassigned"
}
orthogroup_handles = [Patch(color=category_colors[k], label=v) for k, v in orthogroup_labels.items()]

sv_labels = {
    "inversion": "Inversion",
    "insertion_cv": "Deletion",
    "insertion_mg": "Local Insertion",
    "insertion_mg_foreign": "Foreign Insertion",
    "complex_sv": "Complex SV"
}
sv_handles = [Patch(color=sv_colors[k], label=v) for k, v in sv_labels.items()]

# Top orthogroup legend
fig.legend(handles=orthogroup_handles, loc='upper center', title="Orthogroup Category",ncol=len(orthogroup_handles), frameon=False,fontsize=fs,title_fontsize=fs)
# Right-side SV legend
fig.subplots_adjust(right=0.85)  # Make room on the right
legend_ax = fig.add_axes([0.75, 0.25, 0.12, 0.5])  # [left, bottom, width, height]
legend_ax.axis('off')  # Hide frame
legend_ax.legend(handles=sv_handles, loc='center left',title="SV Gap Type", frameon=False, fontsize=fs,title_fontsize=fs)

plt.subplots_adjust(left=0.06, right=0.85, top=0.95, bottom=0.03)
plt.tight_layout(rect=[0.06, 0, 0.85, 0.95])
plt.savefig("orthogroup_chrombar_with_gaps.png", dpi=600)
plt.close()

# sub_gaps = gap_df[
#     (gap_df["cf_chr"] == chrom) &
#     (gap_df["cf_pos1"] < length) & 
#     (gap_df["cf_pos2"] > 0)
# ]
# print(gap_df["cf_chr"].unique())
# print(chrom_sizes["chrom"].unique())
# print("Unique SV types for", chrom, ":", sub_gaps["event_type"].unique())
