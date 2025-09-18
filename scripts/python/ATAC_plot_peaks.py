#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pyBigWig
import pandas as pd
import numpy as np
import os
import json

# ==============================
# USER CONFIGURATION
# ==============================
working_dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2"
os.chdir(working_dir)

ATAC_BIGWIG = "merged_corrected.bw"
FOOTPRINT_BIGWIG = "Footprints_all_peaks.bw"
PEAKS_BED = "consensus_annotated.bed"
GENES_BED = "genes_function.bed"
FOOTPRINT_BED = "SE_footprint_raw.bed"
INS_BM_FILE = "analysis_10000/ALL_10000_TADs_tad_score.bm"

CHROM = "HiC_scaffold_1"
PADDING = 10000
REGION_START = 13283141- PADDING
REGION_END = 13446823 + PADDING

OUTPUT_PNG = f"Regulatory_tracks_{CHROM}_{REGION_START}_{REGION_END}.png"

# ==============================
# LOAD BIGWIG SIGNAL
# ==============================
print("[INFO] Reading ATAC and Footprint bigWigs...")
bw_atac = pyBigWig.open(ATAC_BIGWIG)
bw_foot = pyBigWig.open(FOOTPRINT_BIGWIG)

BIN_SIZE = 200
n_bins = (REGION_END - REGION_START) // BIN_SIZE
bin_edges = np.linspace(REGION_START, REGION_END, n_bins + 1)
bin_centers_mb = (bin_edges[:-1] + bin_edges[1:]) / 2 / 1e6

atac_binned = bw_atac.stats(CHROM, REGION_START, REGION_END, nBins=n_bins)
atac_binned = np.array([0 if v is None else v for v in atac_binned])
bin_positions = np.linspace(REGION_START, REGION_END, n_bins) / 1e6

positions = np.arange(REGION_START, REGION_END)
atac_signal = np.nan_to_num(np.array(bw_atac.values(CHROM, REGION_START, REGION_END)))
foot_signal = np.nan_to_num(np.array(bw_foot.values(CHROM, REGION_START, REGION_END)))

# ==============================
# LOAD BED FILES
# ==============================
print("[INFO] Reading Peaks BED...")
peaks = pd.read_csv(PEAKS_BED, sep="\t", header=None, usecols=[0,1,2,5],
                     names=["chr","start","end","istop"])
peaks = peaks[(peaks["chr"] == CHROM) &
               (peaks["end"] >= REGION_START) &
               (peaks["start"] <= REGION_END)]

print("[INFO] Reading Genes BED...")
genes = pd.read_csv(GENES_BED, sep="\t", header=None, usecols=[0,1,2,3,5],
                     names=["chr","start","end","name","strand"])
genes = genes[(genes["chr"] == CHROM) &
               (genes["end"] >= REGION_START) &
               (genes["start"] <= REGION_END)]

print("[INFO] Reading Footprint BED...")
footprints = pd.read_csv(FOOTPRINT_BED, sep='\t', header=0)

# ==============================
# LOAD INSULATION SCORE
# ==============================
print("[INFO] Reading Insulation BM file...")
with open(INS_BM_FILE) as f:
    header_info = json.loads(f.readline().strip("#").strip())
bm_df = pd.read_csv(INS_BM_FILE, sep='\t', skiprows=1, header=None)
num_scores = bm_df.shape[1] - 3
score_cols = [f"score_{i+1}" for i in range(num_scores)]
bm_df.columns = ['chrom', 'start', 'end'] + score_cols

ins_region = bm_df[(bm_df['chrom'] == CHROM) &
                    (bm_df['end'] >= REGION_START) &
                    (bm_df['start'] <= REGION_END)].copy()
ins_region['mid'] = (ins_region['start'] + ins_region['end']) / 2

region_start_mb = REGION_START / 1e6
region_end_mb = REGION_END / 1e6

# ==============================
# AGGREGATE FOOTPRINTS BY TF CLASS
# ==============================
print("[INFO] Aggregating TF-class binned footprint signal...")
region_footprints = footprints[
    (footprints["chr"] == CHROM) &
    (footprints["end"] >= REGION_START) &
    (footprints["start"] <= REGION_END)
].copy()

region_footprints["TF_class"] = region_footprints["TF_class"].fillna("Unknown")
region_footprints["bin"] = np.searchsorted(bin_edges, region_footprints["start"], side='right') - 1
region_footprints = region_footprints[(region_footprints["bin"] >= 0) & (region_footprints["bin"] < n_bins)]

agg_matrix = (
    region_footprints
    .groupby(["TF_class", "bin"])["footprint_score"]
    .mean()
    .reset_index()
    .pivot(index="bin", columns="TF_class", values="footprint_score")
    .fillna(0)
)

class_sums = agg_matrix.sum(axis=0).sort_values(ascending=False)
MAX_TOP_CLASSES = 4
top_classes = class_sums.head(MAX_TOP_CLASSES).index.tolist()
n_footprint_tracks = len(top_classes)

print(f"[INFO] Top TF classes in region: {top_classes}")
print(f"[INFO] Planning {n_footprint_tracks} dynamic footprint tracks.")

# ==============================
# CREATE PLOT
# ==============================
total_tracks = 1 + n_footprint_tracks + 3
height_ratios = [1.5] + [1.0]*n_footprint_tracks + [0.5, 0.5, 1.0]
fig, axes = plt.subplots(
    nrows=total_tracks,
    ncols=1,
    figsize=(12, max(8, total_tracks * 1.8)),
    sharex=True,
    gridspec_kw={'height_ratios': height_ratios}
)

# ------------------------------
# Track 1: ATAC Signal
# ------------------------------
axes[0].bar(bin_positions, atac_binned, width=BIN_SIZE/1e6, color="black")
axes[0].axhline(0, color='black', linestyle='-', linewidth=0.8)
axes[0].set_ylabel("ATAC signal", fontsize=10)
axes[0].set_title(f"{CHROM}: {REGION_START:,}-{REGION_END:,}", fontsize=12)
axes[0].spines[['top','bottom','right']].set_visible(False)
axes[0].tick_params(axis='x', which='both', bottom=False, labelbottom=False)

# ------------------------------
# Dynamic Footprint Tracks
# ------------------------------
for i, tf_class in enumerate(top_classes):
    ax = axes[1 + i]
    y_signal = np.zeros_like(bin_centers_mb)
    if tf_class in agg_matrix.columns:
        y_signal[agg_matrix.index] = agg_matrix[tf_class].values
    ax.plot(bin_centers_mb, y_signal, color='red')
    ax.set_ylabel(tf_class, fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

# ------------------------------
# Consensus Peaks Track
# ------------------------------
peaks_ax = axes[1 + n_footprint_tracks]
for _, row in peaks.iterrows():
    peak_start_mb = max(row['start'], REGION_START) / 1e6
    peak_end_mb = min(row['end'], REGION_END) / 1e6
    color = 'orange' if row.get('istop', 0) == 1 else 'blue'
    peaks_ax.add_patch(patches.Rectangle(
        (peak_start_mb, 0), peak_end_mb - peak_start_mb, 0.8,
        facecolor=color, alpha=0.5, edgecolor=None
    ))
peaks_ax.set_ylim(0, 1)
peaks_ax.set_ylabel("Consensus\npeaks", fontsize=10)
peaks_ax.spines[['top','right']].set_visible(False)
peaks_ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
peaks_ax.set_yticks([])

# ------------------------------
# Gene Annotation Track
# ------------------------------
genes_ax = axes[2 + n_footprint_tracks]
for _, row in genes.iterrows():
    gene_start_mb = max(row['start'], REGION_START) / 1e6
    gene_end_mb = min(row['end'], REGION_END) / 1e6
    gene_width_mb = gene_end_mb - gene_start_mb
    if gene_width_mb <= 0:
        continue
    y_center = 0.5
    body_height = 0.5
    genes_ax.add_patch(
        patches.Rectangle(
            (gene_start_mb, y_center - body_height/2),
            gene_width_mb, body_height,
            facecolor='green', edgecolor='black', alpha=0.6
        )
    )
    # Arrow head
    arrow_len = 0.002 * (REGION_END - REGION_START) / 1e6
    head_width = 0.12
    if row['strand'] == '+':
        arrow_x = gene_start_mb + (gene_end_mb - gene_start_mb)/2 - arrow_len/2
        genes_ax.arrow(arrow_x, y_center, arrow_len, 0,
                        head_width=head_width, head_length=arrow_len,
                        fc='darkgreen', ec=None, length_includes_head=True)
    else:
        arrow_x = gene_end_mb + (gene_start_mb - gene_end_mb)/2 + arrow_len/2
        genes_ax.arrow(arrow_x, y_center, -arrow_len, 0,
                        head_width=head_width, head_length=arrow_len,
                        fc='darkgreen', ec=None, length_includes_head=True)
    genes_ax.text(
        (gene_start_mb + gene_end_mb) / 2, y_center + body_height/2 + 0.05,
        row['name'], ha='center', va='bottom', fontsize=8
    )
genes_ax.set_ylim(0, 1)
genes_ax.set_ylabel("Genes", fontsize=10)
genes_ax.spines[['top','right']].set_visible(False)
genes_ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
genes_ax.set_yticks([])

# ------------------------------
# Insulation Score Track
# ------------------------------
ins_ax = axes[3 + n_footprint_tracks]
for col in score_cols:
    ins_ax.plot(ins_region['mid'] / 1e6, ins_region[col], lw=1.5, label=col)
ins_ax.set_ylabel("TAD insulation\nscore", fontsize=10)
ins_ax.spines[['top','right']].set_visible(False)
# Set x-axis in Mb with 2 decimal places
xticks = ins_ax.get_xticks()
ins_ax.set_xticks(xticks)  # Explicitly set the tick positions
ins_ax.set_xticklabels([f"{x:.2f}" for x in xticks], fontsize=10)  # Format labels as Mb
ins_ax.set_xlabel("Genomic position (Mb)", fontsize=10)
ins_ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
ins_ax.grid(False)

# ------------------------------
# X-axis
# ------------------------------
axes[-1].set_xlabel("Genomic Position (Mb)", fontsize=10)
axes[-1].set_xlim(region_start_mb, region_end_mb)

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=300)
plt.close()
print(f"[DONE] Plot saved as {OUTPUT_PNG}")