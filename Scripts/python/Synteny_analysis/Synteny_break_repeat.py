#!/usr/bin/env python3

import os
import pandas as pd
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# ================================
# CONFIG
# ================================
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final")

FLANK_SIZE = 50000
INPUT_GAP_CSV       = "gap_records_analysis.csv"
OUTPUT_FLANK_BED    = "breakpoints_cf_inside_gap_windows.bed"
REPEATS_BED         = "all_repeats.bed"
OVERLAP_OUTPUT_BED  = "breakpoint_inside_gap_repeat_overlaps.bed"
CV_GENOME_INDEX     = "primary_dedup_chr_masked_hp_sealed.fa.fai"
BACKGROUND_TILES_BED= f"cv_{FLANK_SIZE}_tiles.bed"

FOCUS_CLASSES = ["LINE", "SINE", "LTR","DNA"]


chrom_map = {f"cf{i}": f"HiC_scaffold_{i}" for i in range(1, 20)}

# ================================
# 1. Load MCScanX gap records
# ================================
df_gaps = pd.read_csv(INPUT_GAP_CSV)
boundary_rows = df_gaps[
    (df_gaps["boundary_label"] == "boundary") & (df_gaps["event_type"] != "trivial")
]

print(f"Loaded {len(df_gaps)} total rows")
print(f"Filtered {len(boundary_rows)} boundary rows with non-trivial events")

# ================================
# 2. Define 100kb *inside-gap* flank windows
# ================================
flank_records = []

for idx, row in boundary_rows.iterrows():
    cf_chr_raw = row["cf_chr"]
    pos1, pos2 = int(row["cf_pos1"]), int(row["cf_pos2"])
    event_type = row["event_type"]
    cf_chr = chrom_map.get(cf_chr_raw, cf_chr_raw)

    break1, break2 = min(pos1, pos2), max(pos1, pos2)

    # 100kb downstream of break1
    start1, end1 = break1, break1 + FLANK_SIZE
    flank_records.append([cf_chr, start1, end1, f"bp_{idx}_downstream_break1", 0, event_type])

    # 100kb upstream of break2
    start2, end2 = max(0, break2 - FLANK_SIZE), break2
    flank_records.append([cf_chr, start2, end2, f"bp_{idx}_upstream_break2", 0, event_type])

print(f"Total directional inside-gap windows defined: {len(flank_records)}")

pd.DataFrame(
    flank_records,
    columns=["chrom", "start", "end", "name", "score", "event_type"]
).to_csv(OUTPUT_FLANK_BED, sep="\t", header=False, index=False)
print(f"*** Saved {OUTPUT_FLANK_BED} ***")

# ================================
# 3. Intersect flanks with repeats
# ================================
breakpoints_bt = pybedtools.BedTool(OUTPUT_FLANK_BED)
repeats_bt = pybedtools.BedTool(REPEATS_BED)

overlaps = breakpoints_bt.intersect(repeats_bt, wa=True, wb=True)
overlaps.saveas(OVERLAP_OUTPUT_BED)
print(f"*** DONE: Saved repeat overlaps to {OVERLAP_OUTPUT_BED} ***")

# ================================
# 4. Parse overlaps and classify repeats
# ================================
overlap_cols = [
    "flank_chrom", "flank_start", "flank_end", "flank_name", "flank_score", "event_type",
    "rep_chrom", "rep_start", "rep_end", "rep_id", "rep_dot", "rep_strand", "rep_classfam"
]

# Ensure all lines have at least 13 fields before creating DataFrame
with open(OVERLAP_OUTPUT_BED) as f:
    lines = [line.strip().split("\t") for line in f if len(line.strip().split("\t")) >= 13]

overlap_df = pd.DataFrame(lines, columns=overlap_cols)

# Convert positions to numeric
for col in ["flank_start", "flank_end", "rep_start", "rep_end"]:
    overlap_df[col] = overlap_df[col].astype(int)

# Clip overlap to flank window
overlap_df["clip_start"] = overlap_df[["flank_start", "rep_start"]].max(axis=1)
overlap_df["clip_end"] = overlap_df[["flank_end", "rep_end"]].min(axis=1)
overlap_df = overlap_df[overlap_df["clip_start"] < overlap_df["clip_end"]].copy()
overlap_df["clip_len"] = overlap_df["clip_end"] - overlap_df["clip_start"]

# Parse top-level repeat class robustly
overlap_df["rep_classfam"] = overlap_df["rep_classfam"].astype(str).str.strip()
overlap_df["rep_class"] = overlap_df["rep_classfam"].str.split("/").str[0].str.strip()

# Focus classes
FOCUS_CLASSES = ["LINE", "SINE", "LTR", "DNA"]
overlap_df = overlap_df[overlap_df["rep_class"].isin(FOCUS_CLASSES)].copy()

print(f"Filtered to {len(overlap_df)} clipped intervals in focus classes")

# ================================
# 5. Summarize non-overlapping bp per flank + class
# ================================
results = []
for (flank_id, event_type, rep_class), group in overlap_df.groupby(["flank_name", "event_type", "rep_class"]):
    intervals = [
        pybedtools.Interval(row["flank_chrom"], int(row["clip_start"]), int(row["clip_end"]))
        for _, row in group.iterrows()
    ]
    merged = pybedtools.BedTool(intervals).sort().merge()
    total_bp = sum(int(x.end) - int(x.start) for x in merged)
    results.append([flank_id, event_type, rep_class, total_bp])

coverage_df = pd.DataFrame(results, columns=["flank_id", "event_type", "rep_class", "repeat_covered_bp"])

# Add zeros for missing combinations
all_flanks = coverage_df[["flank_id", "event_type"]].drop_duplicates()
all_combos = all_flanks.assign(key=1).merge(
    pd.DataFrame({"rep_class": FOCUS_CLASSES, "key":1}), on="key"
).drop("key", axis=1)

coverage_df = all_combos.merge(
    coverage_df, on=["flank_id", "event_type", "rep_class"], how="left"
)
coverage_df["repeat_covered_bp"] = coverage_df["repeat_covered_bp"].fillna(0)

coverage_df.to_csv("foreground_repeat_class_coverage.csv", index=False)
print("*** Saved foreground_repeat_class_coverage.csv ***")

print(list(coverage_df))
print(coverage_df)
# ================================
# 6. Plot: Repeat Class vs Rearrangement Type
# ================================

# ================================
# All types
plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=coverage_df,
    x="rep_class",
    y="repeat_covered_bp",
    hue="event_type",
    showfliers=False,
    order=FOCUS_CLASSES,
    hue_order=["inversion","insertion_cv","insertion_mg","insertion_mg_foreign", "complex_sv" ]  # <-- set your desired order
)

plt.ylabel(f"Repeat-covered bp per {int(FLANK_SIZE/1000)}kb flanking window",fontsize=14)
plt.xlabel("Repeat class",fontsize=14)
plt.ylim(0, 15000)

# Customize legend labels
handles, labels = ax.get_legend_handles_labels()
# Example: Replace labels with custom ones
custom_labels = ['Inversion','Deletion','Local Insertion','Foreign Insertion', 'Complex SV']
ax.legend(handles=handles, labels=custom_labels, title=None, loc='upper right')

plt.tight_layout()
plt.savefig("repeat_class_coverage_by_event_type.png",dpi=600)
plt.close()
print("*** Saved repeat_class_coverage_by_event_type.png ***")

# ================================
# 1 vs rest
from scipy.stats import mannwhitneyu

# Define your direction for each comparison (customize as needed)
comparison_directions = {
    "LINE": "greater",
    "SINE": "greater",
    "LTR": "less",
    "DNA": "greater"
}

# Your custom labels (must match the number/order of unique rep_class in the plot)
custom_labels = ["LINE","SINE","LTR","DNA transposon"]  # example

# Loop through each repeat class
for rep_class, group_df in coverage_df.groupby("rep_class"):
    # Extract insertion_mg group
    group_base = group_df[group_df["event_type"] == "insertion_mg"]["repeat_covered_bp"]
    # Extract all other SV types
    group_compare = group_df[group_df["event_type"] != "insertion_mg"]["repeat_covered_bp"]

    # Get direction for this rep_class
    direction = comparison_directions.get(rep_class, "two-sided")

    if len(group_base) >= 3 and len(group_compare) >= 3:
        stat, p = mannwhitneyu(group_base, group_compare, alternative=direction)
        print(f"{rep_class} | insertion_mg vs Others | alt: {direction:7s} | U = {stat:.2f} | p = {p:.4e}")
    else:
        print(f"{rep_class} | Not enough data (n={len(group_base)} vs {len(group_compare)})")

# Step 1: Create a new column for collapsed event types
coverage_df["event_group"] = coverage_df["event_type"].apply(
    lambda x: "insertion_mg" if x == "insertion_mg" else "Other"
)

# Step 2: Plot the boxplot
plt.figure(figsize=(10, 6))

ax = sns.boxplot(
    data=coverage_df,
    x="rep_class",
    y="repeat_covered_bp",
    hue="event_group",
    showfliers=False,
    order=FOCUS_CLASSES,
    width=0.6,
    hue_order=["insertion_mg", "Other"],  # Order in legend and plot
    # palette={"insertion_mg": "#E41A1C", "Other": "#377EB8"}  # Optional: custom colors
)
# Set y-axis range and limit
y_max = coverage_df["repeat_covered_bp"].max()
y_range = y_max * 1.3
# ax.set_ylim(0, y_range)
plt.ylim(0,15000)
# Bracket annotations
for i, rep_class in enumerate(FOCUS_CLASSES):
    sub_df = coverage_df[coverage_df["rep_class"] == rep_class]
    group1 = sub_df[sub_df["event_group"] == "insertion_mg"]["repeat_covered_bp"]
    group2 = sub_df[sub_df["event_group"] == "Other"]["repeat_covered_bp"]
    direction = comparison_directions.get(rep_class, "two-sided")

    if len(group1) >= 3 and len(group2) >= 3:
        stat, p = mannwhitneyu(group1, group2, alternative=direction)

        # Use 90th percentile as vertical base
        q75_1 = group1.quantile(0.90)
        q75_2 = group2.quantile(0.90)
        this_base = max(q75_1, q75_2)+4000

        # Set bracket height
        bracket_height = this_base + 0.01 * y_range
        line_height = 0.01 * y_range

        # Draw bracket
        ax.plot([i - 0.2, i + 0.2], [bracket_height, bracket_height], color='black', linewidth=1)
        ax.plot([i - 0.2, i - 0.2], [bracket_height - line_height, bracket_height], color='black', linewidth=1)
        ax.plot([i + 0.2, i + 0.2], [bracket_height - line_height, bracket_height], color='black', linewidth=1)

        # Significance label
        if p < 0.001:
            label = "***"
        elif p < 0.01:
            label = "**"
        elif p < 0.05:
            label = "*"
        else:
            label = "ns"
        ax.text(i, bracket_height + 0.01 * y_range, label, ha='center', va='bottom', fontsize=11)

# Step 3: Aesthetics
plt.ylabel(f"Repeat-covered bp per {int(FLANK_SIZE/1000)}kb flanking window", fontsize=14)
plt.xlabel("Repeat class", fontsize=14)  # Give ample space for annotations
# Step 1: get the handles (colors) and labels (keys)
handles, labels = ax.get_legend_handles_labels()

# Step 2: build your own mapping
label_map = {
    "insertion_mg": "Local Insertion",
    "Other": "Other SVs"
}

# Step 3: apply new labels with correct color handles
ax.legend(
    handles=handles,
    labels=[label_map[l] for l in labels],
    title="SV Gap Type",
    loc='upper right'
)

ax.set_xticklabels(custom_labels,fontsize=11)
# plt.subplots_adjust(top=0.92)  # or top=0.9 if you want more space
plt.tight_layout()
plt.savefig("repeat_class_coverage_insertionmg_vs_others.png", dpi=600)
plt.close()
print("*** Saved repeat_class_coverage_insertionmg_vs_others.png ***")



# # ================================
# # 7. Create genome background 100kb tiles
# # ================================
# genome_sizes = pd.read_csv(CV_GENOME_INDEX, sep="\t", header=None, usecols=[0,1], names=["chrom","length"])

# tiles = []
# for _, row in genome_sizes.iterrows():
#     for start in range(0, row["length"], FLANK_SIZE):
#         end = min(start + FLANK_SIZE, row["length"])
#         tiles.append([row["chrom"], start, end])

# pd.DataFrame(tiles, columns=["chrom","start","end"]).to_csv(BACKGROUND_TILES_BED, sep="\t", header=False, index=False)
# print(f"*** Saved {BACKGROUND_TILES_BED} with {len(tiles)} windows ***")

# # ================================
# # 8. Overlap background with repeats
# # ================================
# bg_tiles_bt = pybedtools.BedTool(BACKGROUND_TILES_BED)
# bg_overlaps = bg_tiles_bt.intersect(repeats_bt, wa=True, wb=True)
# bg_overlaps.saveas("background_repeat_overlaps.bed")
# print("*** DONE: Saved background_repeat_overlaps.bed ***")

# # ================================
# # 9. Summarize background coverage per repeat class
# # ================================
# bg_records = []
# for interval in bg_overlaps:
#     fields = interval.fields

#     if len(fields) < 10:
#         continue  # skip malformed lines

#     tile_chrom, tile_start, tile_end = fields[0], fields[1], fields[2]
#     rep_start, rep_end = fields[4], fields[5]
#     rep_classfam = fields[9]

#     # Skip if class/fam missing
#     if rep_classfam == '.' or pd.isna(rep_classfam):
#         continue

#     try:
#         rep_start = int(rep_start)
#         rep_end = int(rep_end)
#         tile_start = int(tile_start)
#         tile_end = int(tile_end)
#     except ValueError:
#         continue

#     rep_class = rep_classfam.split('/')[0].strip()
#     if rep_class not in FOCUS_CLASSES:
#         continue

#     clip_start = max(tile_start, rep_start)
#     clip_end = min(tile_end, rep_end)

#     if clip_start < clip_end:
#         bg_records.append([tile_chrom, tile_start, tile_end, clip_start, clip_end, rep_class])

# bg_df = pd.DataFrame(bg_records, columns=[
#     "tile_chrom","tile_start","tile_end","clip_start","clip_end","rep_class"
# ])
# print(f"Collected {len(bg_df)} clipped background repeat-overlap records")

# # Merge and sum non-overlapping coverage per tile per class
# bg_summaries = []
# for (chrom, start, end, rep_class), group in bg_df.groupby(["tile_chrom","tile_start","tile_end","rep_class"]):
#     intervals = [
#         pybedtools.Interval(chrom, int(row["clip_start"]), int(row["clip_end"]))
#         for _, row in group.iterrows()
#     ]
#     merged = pybedtools.BedTool(intervals).sort().merge()
#     total_bp = sum(int(x.end) - int(x.start) for x in merged)
#     bg_summaries.append([chrom, start, end, rep_class, total_bp])

# bg_summary_df = pd.DataFrame(bg_summaries, columns=["id","start","end","rep_class","repeat_covered_bp"])
# bg_summary_df["group"] = "genome_background"
# bg_summary_df.to_csv("background_repeat_coverage_summary_by_class.csv", index=False)
# print("*** Saved background_repeat_coverage_summary_by_class.csv ***")

# # ================================
# # 10. Combine foreground and background
# # ================================
# fg_df = coverage_df.copy()
# fg_df = fg_df.rename(columns={"flank_id":"id"})
# fg_df["group"] = "breakend_flank"
# print(fg_df)
# combined_df = pd.concat([
#     fg_df[["id","rep_class","repeat_covered_bp","group"]],
#     bg_summary_df[["id","rep_class","repeat_covered_bp","group"]]
# ], ignore_index=True)
# combined_df.to_csv("combined_repeat_coverage_by_class.csv", index=False)
# print("*** Saved combined_repeat_coverage_by_class.csv ***")

# ================================
# 11. Plot grouped boxplot with p-values
# ================================
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

combined_df=pd.read_csv("combined_repeat_coverage_by_class.csv",header=0)

plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=combined_df,
    x="rep_class",
    y="repeat_covered_bp",
    hue="group",
    showfliers=False,
    order=FOCUS_CLASSES,
    width=0.6,
    # palette={"breakend_flank": "#E41A1C", "genome_background": "#377EB8"}  # Optional: custom colors

)
plt.ylabel(f"Repeat-covered bp per {int(FLANK_SIZE/1000)}kb window",fontsize=14)
plt.xlabel("Repeat Class",fontsize=14)
plt.ylim(0,15000)
legend_labels = ['Synteny Boundary', 'Genomic Background']
ax.legend(handles=handles, labels=legend_labels,title="Window Type",loc='upper right')

plt.tight_layout()

# === Compute p-values ===
pval_results = []
for rep_class in FOCUS_CLASSES:
    flank_vals = combined_df[
        (combined_df["group"]=="breakend_flank") & 
        (combined_df["rep_class"]==rep_class)
    ]["repeat_covered_bp"]
    bg_vals = combined_df[
        (combined_df["group"]=="genome_background") & 
        (combined_df["rep_class"]==rep_class)
    ]["repeat_covered_bp"]
    if len(flank_vals) < 2 or len(bg_vals) < 2:
        pval_results.append((rep_class, None))
    else:
        _, p_val = mannwhitneyu(flank_vals, bg_vals, alternative="two-sided")
        pval_results.append((rep_class, p_val))

print("\n=== Mann-Whitney U Tests per Repeat Class ===")
for rep_class, p in pval_results:
    if p is None:
        print(f"{rep_class}: too few values")
    else:
        print(f"{rep_class}: p-value = {p:.4g}")

# === Annotate with brackets ===
ymax = combined_df["repeat_covered_bp"].max()
ymin = combined_df["repeat_covered_bp"].min()
y_range = ymax - ymin

# Get x positions for tick centers
xticks = ax.get_xticks()

for xpos, (rep_class, p_val) in zip(xticks, pval_results):
    # For each group, get data
    subdata = combined_df[combined_df["rep_class"] == rep_class]

    q75_flank = subdata[subdata["group"]=="breakend_flank"]["repeat_covered_bp"].quantile(0.90)
    q75_bg = subdata[subdata["group"]=="genome_background"]["repeat_covered_bp"].quantile(0.90)

    this_base = max(q75_flank, q75_bg)
    bracket_height = this_base + 0.1*y_range
    line_height = 0.01*y_range
    # line_height = 400

    # Plot bracket lines
    ax.plot([xpos - 0.2, xpos + 0.2], [bracket_height, bracket_height], color='black', linewidth=1)
    ax.plot([xpos - 0.2, xpos - 0.2], [bracket_height - line_height, bracket_height], color='black', linewidth=1)
    ax.plot([xpos + 0.2, xpos + 0.2], [bracket_height - line_height, bracket_height], color='black', linewidth=1)

    # Format p-value text
    if p_val is None:
        label = "n/a"
    elif p_val < 0.001:
        label = "***"
    elif p_val < 0.01:
        label = "**"
    elif p_val < 0.05:
        label = "*"
    else:
        label = "ns"

    ax.text(xpos, bracket_height + 0.002*y_range, label, ha='center', va='bottom', fontsize=11)
ax.set_xticklabels(custom_labels,fontsize=11)
plt.savefig("repeat_content_breakend_vs_background_by_class_with_brackets.png",dpi=600)
plt.close()
print("*** Saved repeat_content_breakend_vs_background_by_class_with_brackets.png ***")