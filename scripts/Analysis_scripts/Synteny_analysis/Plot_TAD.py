import os
import cooler
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import json
import pandas as pd
import math

############################################
# === USER SETTINGS ========================
############################################
bin_size = 5000
working_dir = f"/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/analysis_{bin_size}"
os.chdir(working_dir)

cool_file = f'ALL_{bin_size}.cool'

chrom = 'HiC_scaffold_6'
start = 28201094
end = 28321595

region_start = int(start/bin_size)*bin_size
region_end = int(end/bin_size)*bin_size
############################################
# === Load Hi-C Matrix =====================
############################################
print("Loading .cool file...")
clr = cooler.Cooler(cool_file)
print(f"Fetching matrix for region: {chrom}:{region_start}-{region_end}")
matrix = clr.matrix(balance=False).fetch(f'{chrom}:{region_start}-{region_end}')
print(f"Matrix shape: {matrix.shape}")


############################################
# === Log-transform ========================
############################################
matrix_log = np.log1p(matrix)
print(f"Matrix log-transform: min={np.nanmin(matrix_log):.2f}, max={np.nanmax(matrix_log):.2f}")

############################################
# === Rotate matrix by 45° =================
############################################
print("Rotating matrix by 45°...")
rotated=rotate(matrix_log,45,reshape=True,order=1)
nrows,ncols=rotated.shape
center_row=nrows//2

# Only keep the UPPER triangle above diagonal
rotated_upper = rotated[center_row+1:, :]
rotated_upper_masked = np.ma.masked_where(rotated_upper <= 0, rotated_upper)

print(f"Rotated upper shape: {rotated_upper_masked.shape}")

############################################
# === Genomic coordinates for heatmap cols
############################################
num_cols = rotated_upper_masked.shape[1]
bin_edges = np.linspace(start, end, num_cols + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

############################################
# === Load insulation score .bm =============
############################################
bm_file = f'ALL_{bin_size}_TADs_tad_score.bm'

with open(bm_file) as f:
    header_line = f.readline()
header_info = json.loads(header_line.strip("#").strip())

bm_df = pd.read_csv(bm_file, sep='\t', skiprows=1, header=None)
num_scores = bm_df.shape[1] - 3
score_cols = [f"score_{i+1}" for i in range(num_scores)]
bm_df.columns = ['chrom', 'start', 'end'] + score_cols

# Subset to region
ins_region = bm_df[
    (bm_df['chrom'] == chrom) &
    (bm_df['end'] >= start) &
    (bm_df['start'] <= end)
].copy()
ins_region['mid'] = (ins_region['start'] + ins_region['end']) / 2

print(f"Loaded insulation columns: {score_cols}")
print(f"Insulation bins in region: {len(ins_region)}")

############################################
# === Plotting =============================
############################################
print("Plotting...")
fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(12, 8),
    gridspec_kw={'height_ratios': [3, 1]},
    constrained_layout=True
)

# === Top panel: heatmap
cmap = plt.cm.Reds.copy()
cmap.set_bad(color='white')

triangle_height = rotated_upper_masked.shape[0]

adj_fac=2.2
pixel_bp_width = (end - start) / num_cols
heatmap_start = start - adj_fac* pixel_bp_width
heatmap_end = end + adj_fac*pixel_bp_width
extent = [heatmap_start, heatmap_end, 0, triangle_height]

im = ax1.imshow(
    rotated_upper_masked,
    cmap=cmap,
    origin='lower',
    interpolation='none',
    extent=extent,
    aspect='auto'
)

# ax1.set_xlabel(f"{chrom} Position (bp)", fontsize=12)
ax1.yaxis.set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
# ax1.axis('off')

# Add colorbar
cbar = fig.colorbar(im, ax=ax1, fraction=0.03, pad=0.02,shrink=0.3)
cbar.set_label('log(1 + Contact Counts)', fontsize=10)
cbar.ax.tick_params(labelsize=8)



############################################
# === Bottom panel: insulation score track
############################################
for col in score_cols:
    ax2.plot(ins_region['mid'], ins_region[col], lw=1.5, label=col)

# ax2.set_title('Insulation Scores (multiple window sizes)', fontsize=14)
ax2.set_xlabel(f"{chrom} Position (bp)", fontsize=12)
ax2.set_ylabel('Score', fontsize=12)
ax2.grid(False)
ax2.set_xlim(heatmap_start, heatmap_end)
# ax2.legend(fontsize=8, ncol=2)


ax2.yaxis.set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)

############################################
# === Finish ==============================
############################################
# plt.tight_layout()
plt.savefig(f"{bin_size}_{chrom}_{start}_{end}")
plt.show()