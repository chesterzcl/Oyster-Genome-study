import os
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon

# ====== User Inputs ======
fai_file = "primary_dedup_chr_masked_hp_sealed.fa.fai"
gene_file = "gene_density_200k.tsv"
te_file = "repeat_TE_stacked.tsv"
cov_file = "coverage.tsv"

# === Centromere/Telomere Overlay ===
cent_telo_df = pd.read_csv("cent_telo.tsv", sep="\t")
print(cent_telo_df)

# ====== Load Data ======
# Load chromosome lengths from .fai
fai_df = pd.read_csv(fai_file, sep="\t", header=None,usecols=[0, 1], names=["chrom", "length"])
fai_df["chrom"] = fai_df["chrom"].astype(str).str.strip()
chrom_lengths = fai_df.set_index("chrom")["length"].to_dict()

# Only keep chromosomes 1–10
chrom_order = [f"HiC_scaffold_{i}" for i in range(1, 11)]
chrom_lengths = {k: chrom_lengths[k] for k in chrom_order}

# Load gene density and TE content
gene_df = pd.read_csv(gene_file, sep="\t")
te_df = pd.read_csv(te_file, sep="\t")
cov_df = pd.read_csv(cov_file, sep="\t")

# ====== Normalize gene count for color scale ======
# Cap gene_count at 95th percentile for color normalization
cap_value = gene_df["gene_count"].quantile(0.99)
gene_df["gene_count_capped"] = gene_df["gene_count"].clip(upper=cap_value)

# Use the capped values for color mapping
norm_gene = mcolors.Normalize(vmin=gene_df["gene_count_capped"].min(),vmax=gene_df["gene_count_capped"].max())
cmap_gene = plt.cm.get_cmap("coolwarm")


# Coverage
depth_cap = cov_df["depth"].quantile(0.95)
cov_df["depth_capped"] = cov_df["depth"].clip(upper=depth_cap)
norm_cov = mcolors.Normalize(vmin=cov_df["depth_capped"].min(), vmax=cov_df["depth_capped"].max())
cmap_cov = plt.cm.get_cmap("Purples")

# ====== Plotting Setup ======
fig, ax = plt.subplots(figsize=(8, 8))
x_spacing = 2_500_000   # spacing between chromosome columns
bar_width = 1_000_000   # width of each chromosome bar
x_offset = bar_width*1.5
y_offset =0

# ====== Draw Chromosomes One-by-One ======
for chrom in chrom_order:
    chrom_len = chrom_lengths[chrom]

    # Draw chromosome bar
    ax.add_patch(patches.Rectangle((x_offset, y_offset), bar_width, chrom_len,
                                   edgecolor='none', facecolor='whitesmoke',zorder=1))
    ax.add_patch(patches.Rectangle(
    	(x_offset, y_offset),
    	bar_width,
    	chrom_len,
    	facecolor='none',           # no fill
    	edgecolor='black',          # just border
    	linewidth=1.1,
    	zorder=10                   
    	))

    # Gene density inside the chromosome bar
    gene_sub = gene_df[gene_df["chrom"] == chrom]
    for _, row in gene_sub.iterrows():
        color = cmap_gene(norm_gene(row["gene_count"]))
        ax.add_patch(patches.Rectangle(
            (x_offset, row["bin_start"]),
            bar_width,
            row["bin_end"] - row["bin_start"],
            facecolor=color,
            edgecolor=None
        ))


    # ==== COVERAGE (scatter plot, left of chromosome) ====
    cov_sub = cov_df[cov_df["chrom"] == chrom].copy()
    cov_sub["y"] = (cov_sub["start"] + cov_sub["end"]) / 2  # midpoint of bin
    cov_sub["x"] = x_offset - 200_000 - (cov_sub["depth_capped"] / depth_cap) * 500_000

    ax.scatter(cov_sub["x"], cov_sub["y"], 
               color=cmap_cov(norm_cov(cov_sub["depth_capped"])),
               s=10, alpha=0.8, edgecolor='none', label="Coverage" if chrom==chrom_order[0] else None)


    # TE content (right of chromosome)
    te_sub = te_df[te_df["chrom"] == chrom]
    for _, row in te_sub.iterrows():
        bin_y = row["bin_start"]
        bin_h = row["bin_end"] - row["bin_start"]

        scaling_fact=1.4
        retro_w = scaling_fact*bar_width * (row["retro_pct"] / 100)
        dna_w = scaling_fact*bar_width * (row["dna_pct"] / 100)

        # Retrotransposon
        ax.add_patch(patches.Rectangle(
            (x_offset + bar_width + 100_000, bin_y),
            retro_w,
            bin_h,
            facecolor="orange"
        ))

        # DNA transposon
        ax.add_patch(patches.Rectangle(
            (x_offset + bar_width + 100_000 + retro_w + 50_000, bin_y),
            dna_w,
            bin_h,
            facecolor="cyan"
        ))
    # Subset for this chromosome
    markings = cent_telo_df[cent_telo_df["chrom"] == chrom]
    for _, row in markings.iterrows():
        mid_y = (row["start"] + row["end"]) / 2
    # Triangle size and position (to the right of chromosome bar)
        tri_width = bar_width * 0.4
        tri_height = 1_000_000
        x_base = x_offset + scaling_fact*bar_width + 400_000
        if row["type"] == "centromere":
            # Down-pointing triangle
            triangle = Polygon([
            (x_base, mid_y),
            (x_base + tri_width, mid_y - tri_height / 2),
            (x_base + tri_width, mid_y + tri_height / 2)], closed=True, color='darkgreen')
        elif row["type"] == "telomere":
        # Up-pointing triangl
    # Coordinates for a left-pointing triangle
            triangle = Polygon([
            (x_base, mid_y),                                # tip (points left)
            (x_base + tri_width, mid_y - tri_height / 2),   # bottom right
            (x_base + tri_width, mid_y + tri_height / 2)    # top right
            ], closed=True, color='red', zorder=10)
        ax.add_patch(triangle)

    # Add chromosome label
    ax.text(x_offset + bar_width / 2, chrom_len + 2_000_000,
            chrom.replace("HiC_scaffold_", "Chr "), ha='center', fontsize=9)

    x_offset += bar_width + x_spacing

    # # Bar center x
    # bar_center_x = x_offset + bar_width / 2
    # bar_radius_x = bar_width / 2
    # cap_height = bar_width  # make it symmetric (you can reduce to 0.8× if too round)
    # # Top cap (at chrom_len)
    # ax.add_patch(patches.Ellipse(
    # 	(bar_center_x, chrom_len),
    # 	width=bar_width,
    # 	height=cap_height,
    # 	facecolor='whitesmoke',
    # 	edgecolor='black',
    # 	zorder=1
    # 	))

    # # Bottom cap (at y=0)
    # ax.add_patch(patches.Ellipse(
    #     (bar_center_x, 0),
    #     width=bar_width,
    #     height=cap_height,
    #     facecolor='whitesmoke',
    #     edgecolor='black',
    #     zorder=1
    #     ))



# ====== Final Plot Settings ======
ax.set_xlim(-500_000, x_offset)
ax.set_ylim(-500000, max(chrom_lengths.values()) + 5_000_000)
import matplotlib.ticker as ticker

# Y-axis in Mb
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.0f}"))
ax.yaxis.set_major_locator(ticker.MultipleLocator(10_000_000))  # optional
ax.set_ylabel("Genomic Coordinate (Mb)")
ax.set_xticks([])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Gene count colorbar (rightmost)
sm_gene = plt.cm.ScalarMappable(cmap=cmap_gene, norm=norm_gene)
cbar_gene = fig.colorbar(sm_gene, ax=ax, orientation="vertical",shrink=0.3, fraction=0.02, pad=0.01)
cbar_gene.set_label("Gene Density (per 200kb)")

# Coverage depth colorbar (further left)
sm_cov = plt.cm.ScalarMappable(cmap=cmap_cov, norm=norm_cov)
cbar_cov = fig.colorbar(sm_cov, ax=ax, orientation="horizontal",shrink=0.3, fraction=0.02, pad=0.02)
cbar_cov.set_label("WGS Coverage (per 200kb)")

import matplotlib.patches as mpatches

# Define legend patches
retro_patch = mpatches.Patch(color='orange', label='Retrotransposon')
dna_patch = mpatches.Patch(color='cyan', label='DNA Transposon')

# Triangle color legend handles
centro_patch = mpatches.Patch(color='darkgreen', label='Centromere')
telo_patch = mpatches.Patch(color='red', label='Telomere')

# Add to existing legend, or a new one
ax.legend(
    handles=[retro_patch, dna_patch,centro_patch, telo_patch],
    loc='upper right',
    frameon=False,
    title="Repeat Elements",
    fontsize=9,
    title_fontsize=10
)


plt.tight_layout()
plt.savefig("vertical_chromosome_annotation_plot.png", dpi=600)
plt.show()