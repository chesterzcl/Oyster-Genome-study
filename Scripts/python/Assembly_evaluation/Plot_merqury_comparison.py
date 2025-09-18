import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")
# === HiFi-based QV data ===
hifi_data = {
    "Chromosome": [f"Chr{i+1}" for i in range(10)],
    "HiFi_QV": [67.7901, 67.9893, 71.2394, 73.7745, 69.346, 69.1125, 66.5175, 71.4158, 74.1609, 73.3104]
}

# === WGS-based QV data (replace with your actual values) ===
wgs_data = {
    "Chromosome": [f"Chr{i+1}" for i in range(10)],
    "WGS_QV": [43.7, 44.5, 44.4, 43.4, 44.0, 43.8, 44.2, 43.9, 44.5, 43.3]  # <-- Update with real values
}

# Convert to DataFrames
hifi_df = pd.DataFrame(hifi_data)
wgs_df = pd.DataFrame(wgs_data)

# Merge the two DataFrames
df = pd.merge(hifi_df, wgs_df, on="Chromosome")

# Plotting
x = np.arange(len(df["Chromosome"]))  # Label locations
width = 0.35  # Width of bars

fig, ax = plt.subplots(figsize=(8, 6))
bars1 = ax.bar(x - width/2, df["HiFi_QV"], width, label='HiFi QV', color='skyblue')
bars2 = ax.bar(x + width/2, df["WGS_QV"], width, label='WGS QV', color='lightcoral')

# === Add genome-wide average lines ===
hifi_avg = 69.6252
wgs_avg = 43.97
ax.axhline(hifi_avg, color='blue', linestyle='--', linewidth=1.5)
ax.axhline(wgs_avg, color='red', linestyle='--', linewidth=1.5)
ax.text(len(df)-0.5, hifi_avg + 0.5, f'HiFi Avg: {hifi_avg:.2f}', color='blue', fontsize=10)
ax.text(len(df)-0.5, wgs_avg + 0.5, f'WGS Avg: {wgs_avg:.2f}', color='red', fontsize=10)
# Labeling
ax.set_ylabel("QV Score")
ax.set_xlabel("Chromosome")
# ax.set_title("QV per Chromosome: HiFi vs WGS")
ax.set_xticks(x)
ax.set_xticklabels(df["Chromosome"])
ax.legend()
ax.grid(True, linestyle='--', alpha=0.5)

# Add QV values on top of bars (optional)
for bar in bars1 + bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.1f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha='center', va='bottom', fontsize=8)
plt.rcParams["font.family"] = "Arial"
# Remove top and right border (spines)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Turn off grid
ax.grid(False)
plt.tight_layout()
plt.savefig("chromosome_qv_comparison_with_avg.png",dpi=600,bbox_inches='tight')
plt.show()