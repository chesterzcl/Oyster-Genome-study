import pandas as pd
from pycircos import Gcircle, Garc
import os
import collections
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === INPUT FILES ===
fai_file = "primary_dedup_chr_masked_hp_sealed.fa.fai"
coverage_file = "coverage.tsv"
gene_file = "gene_density_200k.tsv"
repeat_file = "repeat_TE_stacked.tsv"
atac_file = "atac_peak_density.tsv"

# === INIT GCIRCLE ===
gc = Gcircle()
fai_df = pd.read_csv(fai_file, sep="\t", header=None,usecols=[0, 1], names=["chrom", "length"])
fai_df["chrom"] = fai_df["chrom"].astype(str).str.strip()
# Add chromosomes
# print(fai_df)
colors = plt.cm.tab10.colors
# colors = ["#f0f0f0", "#e0e0f8", "#f8e0e0", "#e0f8e0", "#f8f0e0","#e0f0f8", "#f8e0f8", "#d0e0e0", "#e8f0e8", "#f0e8f8"]

custom_labels = {
    "HiC_scaffold_1": "Chr1",
    "HiC_scaffold_2": "Chr2",
    "HiC_scaffold_3": "Chr3",
    "HiC_scaffold_4": "Chr4",
    "HiC_scaffold_5": "Chr5",
    "HiC_scaffold_6": "Chr6",
    "HiC_scaffold_7": "Chr7",
    "HiC_scaffold_8": "Chr8",
    "HiC_scaffold_9": "Chr9",
    "HiC_scaffold_10": "Chr10"
}

for idx, row in fai_df.iterrows():
    chrom = row["chrom"]
    length = row["length"]
    
    arc = Garc(
        arc_id=chrom,
        size=row["length"],
        interspace=3,
        raxis_range=(800, 850),
        labelposition=160,
        label_visible=True,
        label=custom_labels.get(chrom, chrom),  # use custom label if available
        facecolor="#f0f0f0",
        edgecolor="black"
    )
    gc.add_garc(arc)


gc.set_garcs(start=2.5,end=357.5)


for _, row in fai_df.iterrows():
    chrom = row["chrom"]
    chrom_len = int(row["length"])


    tick_interval = 10_000_000

    tick_positions = list(range(0, chrom_len + 1, tick_interval))
    tick_labels = [str(int(pos / 1_000_000)) + " Mb" for pos in tick_positions]

    # Tick every 10 Mb, up to chromosome length
    gc.tickplot(
        chrom,
        raxis_range=(850, 860),          # adjust as needed
        tickinterval=10_000_000,
        ticklabels=None               # Show numeric labels
    )

    gc.tickplot(
        chrom,
        tickinterval=10_100_000,
        ticklabels=tick_labels,      # e.g., ["0", "10", "20", ...]
        raxis_range=(820, 830)       # outer ring for labels
    )

    for text in plt.gca().texts:
        if text.get_text().endswith("Mb"):
            text.set_fontsize(5)

    positions = []
    widths = []
    values = []

    for i, start in enumerate(range(0, chrom_len, 5_000_000)):
        end = min(start + 5_000_000, chrom_len)
        if i % 2 == 0:
            # Only include every other 10 Mb bin (e.g., 0–10 Mb, 20–30 Mb...)
            positions.append(start)
            widths.append(end - start)
            values.append(1)

    gc.barplot(
        chrom,
        data=values,
        positions=positions,
        width=widths,
        base_value=0,
        rlim=[0, 1],
        raxis_range=(800, 850),
        facecolor="black",
        spine=False
    )


# === Load GC-content binned data ===
gc_file = "GC.tsv"
gc_df = pd.read_csv(gc_file, sep="\t")
gc_df["width"] = gc_df["end"] - gc_df["start"]
print(gc_df["GC_content"].mean())
print(gc_df["GC_content"].std())

# Define weighted average function
def weighted_avg(group):
    return np.average(group["GC_content"], weights=group["width"])

# Compute weighted average GC per chromosome
gc_per_chrom = gc_df.groupby("chrom").apply(weighted_avg).reset_index(name="weighted_GC")

# Display result
# print(gc_per_chrom)

# === Scale for plotting ===
vmin_gc = gc_df["GC_content"].min()*0.7
vmax_gc = gc_df["GC_content"].max() * 1.05   # add margin

# === Plot per chromosome ===
for chrom, group in gc_df.groupby("chrom"):
    positions = group["start"].tolist()
    widths = group["width"].tolist()
    gc_vals = group["GC_content"].tolist()

    gc.lineplot(
        chrom,
        data=gc_vals,
        positions=positions,
        rlim=[vmin_gc, vmax_gc],
        raxis_range=[740, 790],
        linecolor="black",         # Purple (or use #E41A1C for red)
        linewidth=1.0,
        spine=True
    )


# === 1. Coverage Track (depth) ===
coverage_df = pd.read_csv(coverage_file, sep="\t")

print(coverage_df['depth'].mean())
print(coverage_df['depth'].std())
# === Build arcdata_dict format ===
values_all = []
arcdata_dict = collections.defaultdict(lambda: {"positions": [], "widths": [], "values": []})

for _, row in coverage_df.iterrows():
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    depth = float(row["depth"])
    width = end - start

    arcdata_dict[chrom]["positions"].append(start)
    arcdata_dict[chrom]["widths"].append(width)
    arcdata_dict[chrom]["values"].append(depth)
    values_all.append(depth)

# === Plot coverage as barplot ===
# vmin, vmax = min(values_all), max(values_all)
vmin, vmax = np.percentile(values_all, [1, 99.9])
mean = np.mean(values_all)
std_dev = np.std(values_all, ddof=1)
vmin=mean-2*std_dev
vmax=mean+2*std_dev
print(vmin,vmax)

point=0
clipped=0
for chrom in arcdata_dict:
    values = np.array(arcdata_dict[chrom]["values"])
    positions = np.array(arcdata_dict[chrom]["positions"])
    
    # Clip values outside rlim
    clipped_values = np.clip(values, vmin, vmax)

    # Optional: color clipped points differently
    is_clipped = (values <= vmin) | (values >= vmax)

    total_points = len(values_all)
    point+=total_points
    # Count the number of clipped points
    num_clipped = np.sum(is_clipped)
    clipped+=num_clipped

    colors = ["orangered" if clip else "grey" for clip in is_clipped]

    # Plot all as one scatterplot (optional: assign marker colors)
    gc.scatterplot(
        chrom,
        data=clipped_values.tolist(),
        positions=positions.tolist(),
        rlim=[vmin - 10, vmax + 10],
        raxis_range=(650, 730),
        facecolor=colors,
        spine=True,
        markersize=2.5
    )   


# Calculate percentage
clipped_percentage = (clipped / point) * 100
print(clipped_percentage)

# === Add barplot for repeats ===
repeat_df = pd.read_csv(repeat_file, sep="\t")

values_all = []
arcdata_dict = collections.defaultdict(lambda: {
    "positions": [],
    "widths": [],
    "DNA": [],
    "stacked_Retro": []  # Retro + DNA, for stacked bar base
})

for _, row in repeat_df.iterrows():
    chrom = row["chrom"]
    start = int(row["bin_start"])
    end = int(row["bin_end"])
    width = end - start
    dna_val = float(row["DNA_transposon"])
    retro_val = float(row["Retrotransposon"])

    arcdata_dict[chrom]["positions"].append(start)
    arcdata_dict[chrom]["widths"].append(width)
    arcdata_dict[chrom]["DNA"].append(dna_val)
    arcdata_dict[chrom]["stacked_Retro"].append(dna_val + retro_val)

    values_all.append(dna_val + retro_val)


# === Define rlim ===
vmin = 0
vmax = np.percentile(values_all,99)
# === Add repeat stacked bars ===
for chrom in arcdata_dict:
    positions = arcdata_dict[chrom]["positions"]
    widths = arcdata_dict[chrom]["widths"]
    dna_vals = arcdata_dict[chrom]["DNA"]
    retro_vals = arcdata_dict[chrom]["stacked_Retro"]

    # Compute the height of retro-only part
    retro_only = [retro - dna for retro, dna in zip(retro_vals, dna_vals)]

    # === Plot DNA transposon bars ===
    gc.barplot(
        chrom,
        data=dna_vals,
        positions=positions,
        width=widths,
        base_value=0.0,
        rlim=[vmin, vmax],
        raxis_range=[550, 640],
        facecolor="#F8766D",  # red-orange
        spine=True
    )

    # === Plot retrotransposon (stacked on top of DNA) ===
    for i in range(len(positions)):
        gc.barplot(
            chrom,
            data=[retro_only[i]],
            positions=[positions[i]],
            width=[widths[i]],
            base_value=dna_vals[i],  # base is the height of the DNA bar
            rlim=[vmin, vmax],
            raxis_range=[550, 640],
            facecolor="#00BFC4",  # teal
            spine=False
        )


# === Centromere/Telomere Overlay ===

cent_telo_df = pd.read_csv("cent_telo.txt", sep="\t", header=0, names=["chrom", "start", "end", "type"])

# Overlay colors
overlay_color_map = {
    "centromere": "darkorchid",  # blue
    "telomere": "#e31a1c"     # red
}
# print(cent_telo_df)
# Set the overlay height to slightly higher than repeat bars
overlay_raxis_range = [620, 640]  # Right on top of repeat track
padding=300000

for _, row in cent_telo_df.iterrows():
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    rtype = row["type"].strip().lower()
    color = overlay_color_map.get(rtype, "gray")
    # print(rtype)
    if rtype=="centromere": 
        gc.barplot(
            chrom,
            data=[1],  # dummy height
            positions=[start-padding],
            width=[end - start+2*padding],
            base_value=0.0,
            rlim=[0, 1],  # dummy height scale
            raxis_range=overlay_raxis_range,
            facecolor=color,
            spine=False
        )

    if rtype=="telomere":
        if start>=1_000_000:
            gc.barplot(
                chrom,
                data=[1],  # dummy height
                positions=[start-2*padding],
                width=[end - start+2*padding],
                base_value=0.0,
                rlim=[0, 1],  # dummy height scale
                raxis_range=overlay_raxis_range,
                facecolor=color,
                spine=False
            )
        else:
            gc.barplot(
                chrom,
                data=[1],  # dummy height
                positions=[start],
                width=[end - start+2*padding],
                base_value=0.0,
                rlim=[0, 1],  # dummy height scale
                raxis_range=overlay_raxis_range,
                facecolor=color,
                spine=False
            )                


# === 2. Gene Density Track (with percentile scaling + outlier markers) ===
gene_df = pd.read_csv(gene_file, sep="\t")

# --- Build per-chromosome arcdata dict ---
gene_arcdata = collections.defaultdict(lambda: {"positions": [], "widths": [], "values": []})
all_gene_vals = []

for _, row in gene_df.iterrows():
    chrom  = row["chrom"]
    start  = int(row["bin_start"])
    end    = int(row["bin_end"])
    width  = end - start
    gcount = float(row["gene_count"])

    gene_arcdata[chrom]["positions"].append(start)
    gene_arcdata[chrom]["widths"].append(width)
    gene_arcdata[chrom]["values"].append(gcount)
    all_gene_vals.append(gcount)

# --- Percentile-based limits ---
p5, p95 = np.percentile(all_gene_vals, [1, 99])
buffer  = 0.05 * (p95 - p5)           # 5 % padding
rlim    = [0, p95*2.0] # for barplot scaling

# --- Plot per chromosome ---
for chrom in gene_arcdata:
    vals  = np.array(gene_arcdata[chrom]["values"])
    pos   = np.array(gene_arcdata[chrom]["positions"])
    width = np.array(gene_arcdata[chrom]["widths"])

    # 1. Clip values for the barplot itself
    clipped = np.clip(vals, rlim[0], rlim[1])

    gc.barplot(
        chrom,
        data       = clipped.tolist(),
        positions  = pos.tolist(),
        width      = width.tolist(),
        base_value = 0.0,
        rlim       = rlim,
        raxis_range=[450, 540],   # gene ring
        facecolor  = "#1f77b4",
        spine=True
    )

    # 2. Out-of-range markers (black dots just outside the ring)
    out_mask = (vals < rlim[0]) | (vals > rlim[1])
    if out_mask.any():
        # Push markers ±5 units beyond rlim so they are visible
        marker_vals = np.where(vals[out_mask] > rlim[1],
                               rlim[1] + buffer * 0.2,
                               rlim[0])
        gc.scatterplot(
            chrom,
            data       = marker_vals.tolist(),
            positions  = pos[out_mask].tolist(),
            rlim       = rlim,               # same scaling, but values already offset
            raxis_range=[450, 540],
            facecolor  = "black",
            spine=False
        )

# Load the SV file
# === Load your SD density data ===
sd_file = "SD_density_by_type_100kb.tsv"
df = pd.read_csv(sd_file, sep="\t")

# === Compute bin widths ===
df["width"] = df["end"] - df["start"]

# === Compute bin widths ===
df["SD_bp"] = 0
df["TD_bp"] = 0
# === Compute total bp for scaling ===
df["total_bp"] = df["TD_bp"] + df["SD_bp"] + df["Interchromosomal_bp"]
values_all = df["total_bp"].tolist()

# === Define rlim ===
vmin = 0
vmax = np.percentile(values_all, 95) 
print(f"Barplot rlim = {vmin:.2f} - {vmax:.2f}")

# === Main Circos plotting ===
for chrom, group in df.groupby("chr"):
    positions = group["start"].tolist()
    widths = group["width"].tolist()
    td_vals = group["TD_bp"].tolist()
    sd_vals = group["SD_bp"].tolist()
    inter_vals = group["Interchromosomal_bp"].tolist()

    # Precompute stacking bases
    sd_base = td_vals
    inter_base = [td + sd for td, sd in zip(td_vals, sd_vals)]

    gc.barplot(
        chrom,
        data=td_vals,
        positions=positions,
        width=widths,
        base_value=0.0,
        rlim=[vmin, vmax],
        raxis_range=[390, 440],
        facecolor="#E41A1C",   # Black
        spine=True
    )

    ### --- SD layer (stacked on TD) ---
    for i in range(len(positions)):
        gc.barplot(
            chrom,
            data=[inter_base[i]],
            positions=[positions[i]],
            width=[widths[i]],
            base_value=sd_base[i],
            rlim=[vmin, vmax],
            raxis_range=[390, 440],
            facecolor="#33A02C",
            spine=False
        )

    ### --- Interchromosomal layer (stacked on TD+SD) ---
    # Inter
    for i in range(len(positions)):
        gc.barplot(
            chrom,
            data=[inter_vals[i]],
            positions=[positions[i]],
            width=[widths[i]],
            base_value=inter_base[i], 
            rlim=[vmin, vmax],
            raxis_range=[390, 430],
            facecolor="#984EA3",    # Blue
            spine=False
        )
#################################### === Adding links ===
# link_df=pd.read_csv("SD_link.txt",sep='\t',header=0)
# # print(link_df)
# # Assuming link_df is your DataFrame with the above columns
# for _, row in link_df.iterrows():
#     if row["alignment_length"]>=20000:
#         gc.chord_plot(
#             (row["ref_chr"], row["ref_start"], row["ref_end"],380),
#             (row["qry_chr"], row["qry_start"], row["qry_end"],380),
#             facecolor="#FF8000", 
#             edgecolor="#FF8000",
#             linewidth=1
#         )


#################################### === PLOT ===
gc.figure

# plt.show()
plt.savefig(f"cv_genome_circos.png", dpi=800, bbox_inches='tight')
# plt.close()



# # === 3. ATAC Peak Coverage Track ===
# atac_df = pd.read_csv(atac_file, sep="\t")
# atac_arcdata = collections.defaultdict(lambda: {"positions": [], "widths": [], "values": []})
# all_atac_vals = []

# for _, row in atac_df.iterrows():
#     chrom = row["chrom"]
#     start = int(row["bin_start"])
#     end   = int(row["bin_end"])
#     width = end - start
#     value = float(row["atac_bp"])

#     atac_arcdata[chrom]["positions"].append(start)
#     atac_arcdata[chrom]["widths"].append(width)
#     atac_arcdata[chrom]["values"].append(value)
#     all_atac_vals.append(value)

# # --- Percentile-based limits ---
# p5, p95 = np.percentile(all_atac_vals, [1, 99])
# buffer = 0.05 * (p95 - p5)
# print(p5,p95)
# rlim = [0, p95 + 3500]

# # --- Plot per chromosome ---
# for chrom in atac_arcdata:
#     vals = np.array(atac_arcdata[chrom]["values"])
#     flipped_vals = -vals  # This inverts the bar direction visually
#     pos  = np.array(atac_arcdata[chrom]["positions"])
#     wid  = np.array(atac_arcdata[chrom]["widths"])

#     # 1. Clip for barplot
#     clipped_vals = np.clip(flipped_vals, -rlim[1], -rlim[0])

#     gc.barplot(
#         chrom,
#         data=clipped_vals.tolist(),
#         positions=pos.tolist(),
#         width=wid.tolist(),
#         base_value=0.0,
#         rlim=[-rlim[1], -rlim[0]],
#         raxis_range=[400,480],  # inside of gene track
#         facecolor="#ff7f0e",     # orange
#         spine=True
#     )

#     # 2. Outlier markers (also flipped)
#     out_mask = (vals < rlim[0]) | (vals > rlim[1])
#     if out_mask.any():
#         marker_vals = np.where(vals[out_mask] > rlim[1],
#                                -rlim[1] - buffer * 0.2,
#                                -rlim[0] + buffer * 0.2)

#         gc.scatterplot(
#             chrom,
#             data=marker_vals.tolist(),
#             positions=pos[out_mask].tolist(),
#             rlim=[-rlim[1], -rlim[0]],
#             raxis_range=[400, 480],
#             facecolor="black",
#             spine=False
#         )


####========================================================================================
# Step 1: Merge dataframes on 'chrom'
merged_df = pd.merge(
    gene_df,
    repeat_df[['chrom', 'bin_start', 'bin_end', 'DNA_transposon', 'Retrotransposon']],
    on=['chrom', 'bin_start', 'bin_end']
)
# Step 2: Calculate total repeat content (DNA + Retrotransposon)
merged_df['repeat_total'] = merged_df['DNA_transposon'] + merged_df['Retrotransposon']
# Step 3: Calculate correlation
correlation = merged_df['gene_count'].corr(merged_df['repeat_total'])
corr_dna = merged_df['gene_count'].corr(merged_df['DNA_transposon'])
corr_retro = merged_df['gene_count'].corr(merged_df['Retrotransposon'])

from scipy.stats import pearsonr

# Pearson correlation for DNA transposon
r_dna, p_dna = pearsonr(merged_df['gene_count'], merged_df['DNA_transposon'])

# Styling
sns.set(style="ticks", context="notebook", font_scale=1.2)

plt.figure(figsize=(6, 5))
ax = sns.regplot(
    data=merged_df,
    x='DNA_transposon',
    y='gene_count',
    scatter_kws={'alpha': 0.5, 's': 30},
    line_kws={'color': 'darkred'}
)

# Axis labels
ax.set_xlabel("DNA Transposon Content (bp)")
ax.set_ylabel("Gene Count")

# Add r and p-value text inside plot
plt.text(
    0.05, 0.95,
    f"r = {r_dna:.2f}\np = {p_dna:.1e}",
    transform=ax.transAxes,
    verticalalignment='top',
    horizontalalignment='left',
    fontsize=11,
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray')
)

# Remove gridlines
ax.grid(False)
plt.tight_layout()
plt.savefig("panel4A_gene_vs_dna_transposon.png", dpi=300)
plt.show()



# Pearson correlation for retrotransposon
r_retro, p_retro = pearsonr(merged_df['gene_count'], merged_df['Retrotransposon'])

plt.figure(figsize=(6, 5))
ax = sns.regplot(
    data=merged_df,
    x='Retrotransposon',
    y='gene_count',
    scatter_kws={'alpha': 0.5, 's': 30},
    line_kws={'color': 'darkred'}
)

# Axis labels
ax.set_xlabel("Retrotransposon Content (bp)")
ax.set_ylabel("Gene Count")

# Add r and p-value inside plot
plt.text(
    0.05, 0.95,
    f"r = {r_retro:.2f}\np = {p_retro:.1e}",
    transform=ax.transAxes,
    verticalalignment='top',
    horizontalalignment='left',
    fontsize=11,
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray')
)

# Remove gridlines
ax.grid(False)
plt.tight_layout()
plt.savefig("panel4B_gene_vs_retrotransposon.png", dpi=300)
plt.show()

####========================================================================================



