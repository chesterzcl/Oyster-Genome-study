import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === Data ===
busco_data = {
    "Category": ["Single-copy", "Duplicated", "Fragmented", "Missing"],
    "C.virginica_3.0": [79.4, 18.7, 0.3, 1.6],
    "New": [98.6, 1.1, 0.1, 0.2]
}
df = pd.DataFrame(busco_data)
df.set_index("Category", inplace=True)
df = df.T  # Each genome becomes a bar

# === Plot settings ===
colors = {
    "Single-copy": "#1b9e77",     # Green
    "Duplicated": "#7570b3",      # Purple
    "Fragmented": "#d95f02",      # Orange
    "Missing": "#e7298a"          # Pink
}
category_order = ["Single-copy", "Duplicated", "Fragmented", "Missing"]

# === Plot ===
fig, ax = plt.subplots(figsize=(8, 6))  # Size for 1-column figure
bottoms = [0] * len(df)

# Store y positions of previous side labels to avoid collisions
side_label_y = {i: [] for i in range(len(df))}

# Before plotting, enforce stacking order
category_order = ["Single-copy", "Duplicated", "Fragmented", "Missing"]
df = df[category_order]  # Reorder columns in desired stacking order

# Then proceed with:
for cat in category_order:
    values = df[cat]
    ax.bar(
        df.index,
        values,
        bottom=bottoms,
        color=colors[cat],
        width=0.4,     # or 0.6 depending on your preference
        label=cat
    )

    for i, (x, y, b) in enumerate(zip(df.index, values, bottoms)):
        if y > 1:
            ax.text(i, b + y / 2, f"{y:.1f}%", ha="center", va="center", fontsize=10)
        elif y > 0:
            # === Side label logic for small values ===
            x_offset = 0.2
            # Key change: base_y is now independent of y, more spaced from the top
            base_y = b + y + 0.5  # slightly above the bar segment
            y_pos = base_y

            # Optional collision check (still safe to include)
            while any(abs(y_pos - yy) < 1.5 for yy in side_label_y[i]):
                y_pos += 1.2

            side_label_y[i].append(y_pos)

            # Tick line
            ax.plot([i + x_offset, i + x_offset + 0.1], [y_pos, y_pos], color="black", lw=0.8)

            # Side label
            ax.text(
                i + x_offset + 0.12,
                y_pos,
                f"{cat}: {y:.1f}%",
                ha="left",
                va="center",
                fontsize=10
            )

    # Update stacking baseline
    bottoms = [sum(x) for x in zip(bottoms, values)]

# === Final touches ===
ax.set_ylabel("BUSCO %", fontsize=11)
ax.set_xlabel("Assembly", fontsize=11)
# ax.set_title("BUSCO Composition of Genome Assemblies", fontsize=12)
ax.set_ylim(0, 105)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df.index, fontsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(
    title="BUSCO Category",
    fontsize=10,
    title_fontsize=11,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),  # right of plot
    frameon=False
)
plt.rcParams["font.family"] = "Arial"
plt.tight_layout()
plt.savefig("busco_stacked_barplot_pubgrade.png", dpi=600)
plt.show()