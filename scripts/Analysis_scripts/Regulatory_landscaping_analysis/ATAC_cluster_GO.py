import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Set working directory if needed
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === Load your GO enrichment data ===
high = pd.read_csv("GO_enrichment_high_var.csv")
low = pd.read_csv("GO_enrichment_low_var.csv")

# === Filter for Biological Process terms ===
high_bp = high[high["source"] == "GO:BP"].copy()
low_bp = low[low["source"] == "GO:BP"].copy()

# === Compute enrichment score and -log10(p-value) ===
for df in [high_bp, low_bp]:
    df["enrichment_score"] = (df["intersection_size"] / df["query_size"]) / \
                             (df["term_size"] / df["effective_domain_size"])
    df["-log10_pval"] = -np.log10(df["p_value"])

# === Get top 10 unique GO terms from both groups ===
top_high = high_bp.nlargest(10, "-log10_pval")
top_low = low_bp.nlargest(10, "-log10_pval")

top_terms = pd.concat([top_high, top_low], ignore_index=True)[["native", "name"]].drop_duplicates()

# === Build merged table manually ===
records = []
for _, row in top_terms.iterrows():
    term_id = row["native"]
    term_name = row["name"]

    # Look up in both datasets
    high_row = high_bp[high_bp["native"] == term_id]
    low_row = low_bp[low_bp["native"] == term_id]

    # Extract info
    if not high_row.empty:
        records.append({
            "name": term_name,
            "group": "High variance",
            "enrichment_score": high_row["enrichment_score"].values[0],
            "-log10_pval": high_row["-log10_pval"].values[0]
        })
    else:
        records.append({
            "name": term_name,
            "group": "High variance",
            "enrichment_score": np.nan,
            "-log10_pval": np.nan
        })

    if not low_row.empty:
        records.append({
            "name": term_name,
            "group": "Low variance",
            "enrichment_score": low_row["enrichment_score"].values[0],
            "-log10_pval": low_row["-log10_pval"].values[0]
        })
    else:
        records.append({
            "name": term_name,
            "group": "Low variance",
            "enrichment_score": np.nan,
            "-log10_pval": np.nan
        })

df = pd.DataFrame(records)

# === Set x/y for plotting ===
df["x"] = df["group"].map({"High variance": -0.25, "Low variance": 0.25})
term_order = pd.Series(df["name"].unique()).sort_values().tolist()
# Sort GO terms by –log10(pval) in High variance group
term_order = (
    df[df["group"] == "High variance"]
    .sort_values("-log10_pval", ascending=True)["name"]
    .drop_duplicates()
    .tolist()
)

# Preserve full term list (in case Low-only terms exist)
remaining_terms = df[~df["name"].isin(term_order)]["name"].unique().tolist()
term_order += [t for t in remaining_terms if t not in term_order]

# Assign y positions
df["y"] = df["name"].map({term: i for i, term in enumerate(term_order)})

# === Plot ===
plt.figure(figsize=(9, 10))
sns.set(style="white")


# Custom mapping function for size (adjust for contrast)
def size_map(val):
    if val <= 10:
        return 3 + (val ** 2.5) * 1.3
    else:
        base = 3 + (10 ** 2.5) * 1.3  # ≈ 161.1
        return base + ((val - 10) ** 2.5) * 0.35

# === Main scatterplot ===
ax = sns.scatterplot(
    data=df,
    x="x",
    y="y",
    size="-log10_pval",
    hue="enrichment_score",
    palette="Blues",
    sizes=(size_map(df["-log10_pval"].min()), size_map(df["-log10_pval"].max())),
    edgecolor=None,
    alpha=1,
    legend=False
)

legend_vals = [1, 5, 10]
handles = [
    plt.scatter([], [], s=size_map(v), color='gray', alpha=0.6, label=str(v))
    for v in legend_vals
]

legend = ax.legend(
    handles=handles,
    title='–log10(p-value)',
    loc='lower center',
    bbox_to_anchor=(0.5, -0.15),
    frameon=False,  # removes the border
    fancybox=False,  # makes sure no rounded border
    edgecolor=None,  # avoids any line
    ncol=3,
    handletextpad=0.8,
    columnspacing=1.0,
    fontsize=10,
    title_fontsize=12
)

# === Continuous colorbar for enrichment score ===
norm = plt.Normalize(df["enrichment_score"].min(), df["enrichment_score"].max())
sm = plt.cm.ScalarMappable(cmap="Blues", norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
cbar.set_label("Fold enrichment", fontsize=12)


# === Axis formatting ===
ax.set_xticks([-0.25, 0.25])
ax.set_xticklabels(["High variance peaks", "Low variance peaks"], fontsize=10)
ax.set_xlim(-0.6, 0.6)
ax.set_yticks(df["y"])
ax.set_yticklabels(df["name"], fontsize=10)
ax.set_xlabel("")
ax.set_ylabel("GO Biological Process")
ax.set_title("GO:BP Enrichment — High vs. Low Variance Peaks", fontsize=14)
ax.grid(False)

plt.savefig("Pathway_enrichment_bubble_plot.png")
plt.tight_layout()
plt.show()