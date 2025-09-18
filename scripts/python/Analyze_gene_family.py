import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import from_indicators, UpSet,from_memberships
import os

# -------------------------------------------------------------------
# 1. LOAD THE ORTHOGROUP GENE-COUNT TABLE
# -------------------------------------------------------------------

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/Orthogroups")
file_path="Orthogroups.GeneCount.annotated.tsv"     
df=pd.read_csv(file_path, sep=r"\t")

# Columns holding species counts
species_cols = ["CV","MA","MG","OE"]

# -------------------------------------------------------------------
# 2. QUICK SUMMARY
# -------------------------------------------------------------------
print(f"Rows (orthogroups): {len(df)}")
print("\nTop 5 orthogroups by total count:\n")
print(df.nlargest(5, "Total"))

# -------------------------------------------------------------------
# 3. BAR-PLOT — TOP 20 BY TOTAL
# -------------------------------------------------------------------
top_n=20
top_df=df.nlargest(top_n, "Total").set_index("Orthogroup")

plt.figure(figsize=(10, 5))
plt.bar(top_df.index, top_df["Total"], color="#2ca02c")
plt.title(f"Top {top_n} Orthogroups by Total Gene Count")
plt.ylabel("Gene count (all species)")
plt.xticks(rotation=70, ha="right")
plt.tight_layout()
plt.savefig("barplot_top20_orthogroups.png", dpi=300)
plt.close()

# -------------------------------------------------------------------
# 4. HEATMAP — TOP 50 ORTHOGROUPS
# -------------------------------------------------------------------
top50=df.nlargest(50,"Total").set_index("Orthogroup")

fig,ax=plt.subplots(figsize=(7, 9))
im=ax.imshow(top50[species_cols],aspect="auto",cmap="YlGn")

cbar=fig.colorbar(im,ax=ax)
cbar.set_label("Gene count",rotation=270,labelpad=15)

ax.set_xticks(range(len(species_cols)))
ax.set_xticklabels(species_cols, rotation=45, ha="right")
ax.set_yticks(range(len(top50)))
ax.set_yticklabels(top50.index)
ax.set_title("Heatmap of Gene Counts (Top 50 Orthogroups)")
plt.tight_layout()
plt.savefig("heatmap_top50_orthogroups.png", dpi=300)
plt.close()

# -------------------------------------------------------------------
# 5. UPSET PLOT — PRESENCE/ABSENCE with CUSTOM ORDER
# -------------------------------------------------------------------

# Convert counts to presence/absence (True/False)
presence = df.copy()
for col in species_cols:
    presence[col] = presence[col] > 0

# Reorder the columns according to desired species order
desired_order = ['MA','MG','CV','OE']
presence = presence[desired_order]

# Create UpSet data from indicators
upset_data = from_indicators(desired_order, presence)

# Plot
plt.figure(figsize=(9, 5))
UpSet(
    upset_data,
    subset_size='count',
    show_counts=True,
    sort_by='degree',           # you can also try 'cardinality'
    sort_categories_by=None,    # disable auto-sorting of species
    orientation="horizontal"
).plot()

plt.tight_layout()
plt.savefig("upset_orthogroups_ordered.svg", dpi=600)
plt.close()