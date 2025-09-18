import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final")

# === 1. Load your table ===
df = pd.read_csv("FRRS1_domains.txt",header=0, sep="\t")  # replace with your file
# Expected columns: Gene_id, Total length, Begin, End, Domain

# Rename the domain for plotting
df['Domain'] = df['Domain'].replace(
    {"Cytochrome b561/ferric reductase transmembrane": "Cytochrome b561 domain"}
)
# === 2. Define custom y-axis order ===
custom_order = ["g488","g954","g952", "g950", "g949"]  # Change to your preferred order
df['Gene_id'] = pd.Categorical(df['Gene_id'], categories=custom_order, ordered=True)
genes = custom_order

# === 3. Domain colors ===
domain_colors = {
    "Reeler domain": "#E41A1C",
    "DOMON domain": "#377EB8",
    "Cytochrome b561 domain": "#4DAF4A"
}

# === 4. Create plot ===
plt.figure(figsize=(12, len(genes) * 0.8))
ax = plt.gca()

for i, gene in enumerate(genes):
    sub = df[df['Gene_id'] == gene].sort_values('Begin')
    if sub.empty:
        continue
    total_len = sub['Total length'].iloc[0]

    # Compute segments of the backbone line between domains
    last_end = 0
    for _, row in sub.iterrows():
        start, end = row['Begin'], row['End']

        # Draw line segment before this domain
        ax.plot([last_end, start], [i, i], color='black', lw=2)
        last_end = end

    # Draw line segment after the last domain
    if last_end < total_len:
        ax.plot([last_end, total_len], [i, i], color='black', lw=2)

    # Draw domains (with no internal line)
    for _, row in sub.iterrows():
        start, end, dom = row['Begin'], row['End'], row['Domain']
        color = domain_colors.get(dom, "#999999")
        rect = patches.Rectangle((start, i - 0.2), end - start, 0.4,
                                 facecolor=color, edgecolor='black')
        ax.add_patch(rect)
        # Add domain label
        ax.text((start + end) / 2, i + 0.3, dom, ha='center', va='bottom', fontsize=8)

# === 5. Styling ===
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=10)
ax.tick_params(axis='y', length=0)  # Hide Y-axis ticks but keep labels

ax.set_xlabel("Protein Length (aa)")
ax.set_ylabel("Gene")
# ax.set_title("Protein Domain Architecture of Gene Family", pad=15)

# Remove top/right borders
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.ylim(-1, len(genes))
plt.tight_layout()
plt.savefig("FRRS1_domains.png",dpi=600)