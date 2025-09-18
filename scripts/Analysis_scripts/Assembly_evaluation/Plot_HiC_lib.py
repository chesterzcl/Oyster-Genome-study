#!/usr/bin/env python

import cooler
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 


########################################
# === USER SETTINGS ====================
########################################
binsize = 1000  # adjust if needed
cool_file = f'ALL_{binsize}.cool'
os.chdir(f"/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/analysis_{binsize}")
min_distance = binsize
max_distance = 20_000_000
num_bins = 100

########################################
# === LOAD .cool FILE ==================
########################################
print(f"Loading .cool file: {cool_file}")
clr = cooler.Cooler(cool_file)
bins = clr.bins()[:]
pixels = clr.pixels()[:]
print(f"Loaded: {len(pixels)} contact pixels")

########################################
# === MERGE BIN COORDINATES ===========
########################################
pixels = pixels.merge(bins.add_prefix('bin1_'), left_on='bin1_id', right_index=True)
pixels = pixels.merge(bins.add_prefix('bin2_'), left_on='bin2_id', right_index=True)

########################################
# === FILTER CIS-CONTACTS ==============
########################################
cis_pixels = pixels[pixels['bin1_chrom'] == pixels['bin2_chrom']].copy()
print(f"Cis contacts: {len(cis_pixels)}")

########################################
# === COMPUTE DISTANCES ================
########################################
cis_pixels['distance'] = np.abs(cis_pixels['bin1_start'] - cis_pixels['bin2_start'])

# Filter out zero distance (self-interactions if undesired)
cis_pixels = cis_pixels[cis_pixels['distance'] > 0]
print(f"Non-zero distance contacts: {len(cis_pixels)}")

########################################
# === BINNING (LOG-SPACED) =============
########################################
distance_bins = np.logspace(
    np.log10(min_distance),
    np.log10(max_distance),
    num_bins
)

hist_counts, bin_edges = np.histogram(
    cis_pixels['distance'],
    bins=distance_bins,
    weights=cis_pixels['count']
)

bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

########################################
# === SAVE BINNED DATA =================
########################################
out_table = pd.DataFrame({
    'distance_bp': bin_centers.astype(int),
    'contact_count': hist_counts.astype(int)
})
out_table.to_csv('distance_decay_counts.tsv', sep='\t', index=False)
print("Saved binned data to: distance_decay_counts.tsv")

########################################
# === PLOTTING ========================
########################################
plt.figure(figsize=(6,4))

plt.loglog(
    bin_centers,
    hist_counts,
    marker='o',
    linestyle='-',
    color='black',
    lw=1.5,
    markersize=5,
    label='Contact Decay'
)

plt.xlabel('Genomic Distance (bp)', fontsize=12)
plt.ylabel('Contact Count', fontsize=12)
plt.title('Hi-C Contact Distance Decay', fontsize=14, pad=10)

plt.grid(True, which="both", linewidth=0.5, alpha=0.7)
plt.legend(fontsize=10, loc='upper right')

# Tidy layout
plt.tight_layout()

# Save figure
plt.savefig('distance_decay_plot.pdf')
plt.savefig('distance_decay_plot.png', dpi=300)
plt.show()

print("Saved plot to: distance_decay_plot.pdf / .png")