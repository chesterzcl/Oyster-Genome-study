import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Load HOMER output ===
COMP_NAME="se_vs_distal"
os.chdir(f"/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/{COMP_NAME}")
df = pd.read_csv("knownResults.txt", sep='\t')

# === Parameters ===
LABEL_MIN_LOG10P = 20           # Minimum -log10(P-value) to label
LABEL_MIN_FE =1.5
LABEL_MIN_LOG2FE = np.log2(LABEL_MIN_FE) # Minimum log2(FoldEnrichment) to label


# === Clean and convert columns ===
df["Log P-value"] = df["Log P-value"].astype(float)
df["% of Target Sequences with Motif"] = df["% of Target Sequences with Motif"].str.rstrip('%').astype(float)
df["% of Background Sequences with Motif"] = df["% of Background Sequences with Motif"].str.rstrip('%').astype(float)

# === Calculate statistics ===
df["FoldEnrichment"] = df["% of Target Sequences with Motif"] / df["% of Background Sequences with Motif"]
df["log2FE"] = np.log2(df["FoldEnrichment"])

# === Convert -ln(P) to -log10(P) ===
df["log10_p"] = -df["Log P-value"] / np.log(10)  # ln(P) → log10(P)

# === Bonferroni correction (approximate threshold on log10 scale) ===
bonferroni_log10_cutoff = -np.log10(0.05 / len(df))
df["Significant"] = df["log10_p"] >= bonferroni_log10_cutoff

# === Extract protein name ===
df["Protein"] = df["Motif Name"].str.extract(r"^([^(/]+)")

# === Label selected points ===
df["Label"] = df.apply(
    lambda row: row["Protein"] if (row["log10_p"] >= LABEL_MIN_LOG10P and row["log2FE"] >= LABEL_MIN_LOG2FE) else "",
    axis=1
)

# === Volcano plot ===
plt.figure(figsize=(6,8))
sns.scatterplot(
    data=df,
    x="FoldEnrichment",
    y="log10_p",
    hue="log10_p",
    palette="Reds",
    edgecolor=None,
    alpha=1,
    s=10
)

top_hits = df[
    (df["FoldEnrichment"] >= LABEL_MIN_FE) &
    (df["log10_p"] >= LABEL_MIN_LOG10P)
].sort_values("log10_p", ascending=False).head(10)
# Add labels
from adjustText import adjust_text

texts = []
for _, row in top_hits.iterrows():
    texts.append(
        plt.text(
            row["FoldEnrichment"],
            row["log10_p"],
            row["Protein"],
            ha='left',
            va='bottom',
            weight='bold',
            fontsize=6
        )
    )

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5),force_text=0.4,only_move={'points':'y', 'text':'y'})
# Axes and title
plt.xlabel("Fold Enrichment",fontsize=13)
plt.ylabel("-log10(P-value)",fontsize=13)
xticks = np.arange(0,3,0.5)  # adjust range as needed
plt.xticks(xticks)
plt.xlim(0,3)
# plt.axhline(y=LABEL_MIN_LOG10P, color='black', linestyle='--', linewidth=1)
plt.axvline(x=LABEL_MIN_FE, color='black', linestyle='--', linewidth=1)
# plt.title("DNA-binding motif enrichment in ATAC-seq peaks")
# plt.legend(title="Significant after Bonferroni correction",loc="upper left")
plt.legend().remove()  # for the current figure
plt.tight_layout()
plt.savefig(f"{COMP_NAME}.png",dpi=300)
plt.show()