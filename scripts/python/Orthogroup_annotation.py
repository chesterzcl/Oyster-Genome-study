import os
import re
from collections import Counter
import pandas as pd
import numpy as np
import seaborn as sns
from statsmodels.stats.multitest import multipletests

# Set working directory
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/Orthogroups/")


from Bio import Phylo
import matplotlib.pyplot as plt

# Load the tree with branch lengths
tree = Phylo.read("SpeciesTree_rooted.txt", "newick")

# Function to label branches safely
def safe_branch_label(clade):
    if clade.branch_length is None:
        return ""
    return f"{clade.branch_length:.3f}"


def parse_cv_gff(gff_file):
    gene_names = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attributes = parts[8]
            gene_id_match = re.search(r'ID=([^;]+)', attributes)
            desc_match = re.search(r'Description=([^;]+)', attributes)
            name_match = re.search(r'Name=([^;]+)', attributes)

            if gene_id_match:
                gene_id = gene_id_match.group(1)
                desc = desc_match.group(1) if desc_match else ""
                name = name_match.group(1) if name_match else ""

                # if desc and "uncharacterized" not in desc.lower():
                #     gene_names[gene_id] = desc
                # elif name:
                #     gene_names[gene_id] = name
                # else:
                #     gene_names[gene_id] = ""
                if name:
                    gene_names[gene_id] = name
                elif desc and "uncharacterized" not in desc.lower():
                    gene_names[gene_id] = desc
                else:
                    gene_names[gene_id] = ""
                # if "uncharacterized" in desc.lower():
                #     gene_names[gene_id] = ""
                # else:
                #     gene_names[gene_id] = desc
    return gene_names

def parse_other_gff(gff_file):
    gene_names = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attributes = parts[8]
            gene_id_match = re.search(r'ID=([^;]+)', attributes)
            desc_match = re.search(r'description=([^;]+)', attributes)
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                desc = desc_match.group(1) if desc_match else ""
                if "uncharacterized" in desc.lower():
                    gene_names[gene_id] = ""
                else:
                    gene_names[gene_id] = desc
    return gene_names

def annotate_orthogroups(orthogroup_file, output_file, mapping_dict):
    og_to_annotation = {}

    with open(orthogroup_file) as infile, open(output_file, "w") as outfile:
        for line in infile:
            if ":" not in line:
                outfile.write(line)
                continue

            group, genes_str = line.strip().split(":", 1)
            genes = genes_str.strip().split()

            annotations = []
            annotated_genes = []
            for g in genes:
                species, gene_id = g.split("|")
                name = mapping_dict.get(gene_id, "")
                if name:
                    annotations.append(name)
                annotated_genes.append(f"{g}({name})" if name else f"{g}( )")

            # Get most frequent non-empty annotation
            if annotations:
                freq_anno = Counter(annotations).most_common(1)[0][0]
                og_to_annotation[group] = freq_anno
            else:
                og_to_annotation[group] = ""

            outfile.write(f"{group}: {' '.join(annotated_genes)}\n")

    return og_to_annotation

# === Run ===

cv_map = parse_cv_gff("CV.gff3")
oe_map = parse_other_gff("OE.gff")
# me_map = parse_other_gff("ME.gff")
ma_map = parse_other_gff("MA.gff")
mg_map = parse_other_gff("MG.gff")

# gene_map = {**cv_map, **me_map, **ma_map, **mg_map, **oe_map}
gene_map = {**cv_map, **ma_map, **mg_map, **oe_map}
# gene_map = {**cv_map,**mg_map}

og_to_annotation = annotate_orthogroups("Orthogroups.txt", "Orthogroups.annotated.txt", gene_map)

# Annotate gene count table
count_df = pd.read_csv("Orthogroups.GeneCount.tsv", sep="\t")
count_df["Annotation"] = count_df["Orthogroup"].map(og_to_annotation)
count_df=count_df.rename(columns={"CV_proteins_clean":"CV","MG_proteins_clean":"MG","MA_proteins_clean":"MA","OE_proteins_clean":"OE"})
count_df.to_csv("Orthogroups.GeneCount.annotated.tsv", sep="\t", index=False)






# ---- 0) Start from your dataframe ----
df = count_df.copy()

# --------------------------
# Clean & preprocess
# --------------------------
# Fill NA and calculate Magallana total
df[["CV", "MA", "MG"]] = df[["CV", "MA", "MG"]].fillna(0).astype(int)
df["Magallana"] = df["MA"] + df["MG"]

# Compute total genes across all orthogroups
total_cv = df["CV"].sum()
total_mgma = df["Magallana"].sum()

# --------------------------
# Compute log2 Fold Change
# --------------------------
df["log2FC"] = np.log2((df["CV"] + 1) / (df["Magallana"] + 1))

# --------------------------
# Fisher's exact test
# --------------------------
def fisher_pval(row):
    # 2x2 contingency table:
    #               In OG     Not in OG
    # CV           a = row["CV"]     b = total_cv - row["CV"]
    # Magallana    c = row["Magallana"]  d = total_mgma - row["Magallana"]
    a = row["CV"]
    b = total_cv - a
    c = row["Magallana"]
    d = total_mgma - c
    table = [[a, b], [c, d]]
    try:
        _, p = fisher_exact(table, alternative="two-sided")
    except:
        p = np.nan
    return p

df["pval"] = df.apply(fisher_pval, axis=1)

# --------------------------
# Adjust p-values (FDR)
# --------------------------
df["fdr"] = multipletests(df["pval"], method="fdr_bh")[1]

# --------------------------
# Save result
# --------------------------
output_cols = ["Orthogroup", "CV", "MA", "MG", "Magallana", "log2FC", "pval", "fdr", "Annotation"]
df[output_cols].to_csv("orthogroups_with_log2FC_fdr.tsv", sep="\t", index=False)

# --------------------------
# Summary: Top expansions and contractions
# --------------------------
expanded_df = df[(df["log2FC"] > 1) & (df["fdr"] < 0.05)]
contracted_df = df[(df["log2FC"] < -1) & (df["fdr"] < 0.05)]

print(f"Significant CV expansions: {len(expanded_df)} orthogroups")
print(f"Significant CV contractions: {len(contracted_df)} orthogroups")

# Optional: Save these too
expanded_df.to_csv("significant_cv_expansions.tsv", sep="\t", index=False)
contracted_df.to_csv("significant_cv_contractions.tsv", sep="\t", index=False)


# Classification
FOLD = 2
def classify_row(r):
    cv, mgma = r["CV"], r["Magallana"]
    if (cv >= FOLD * mgma) and (mgma >= 1):
        return "CV_expanded"
    elif (mgma >= FOLD * cv) and (cv >= 1):
        return "CV_contracted"
    elif (cv > 0) and (mgma == 0):
        return "CV_specific"
    elif (cv == 0) and (mgma > 0):
        return "Magallana_specific"
    elif (cv > 0) and (mgma > 0) and (1/FOLD <= cv / (mgma + 1e-9) <= FOLD):
        return "Conserved_like"
    else:
        return "Other"

df["Category"] = df.apply(classify_row, axis=1)

# Ratio for sorting
eps = 1e-9
df["cv_vs_mgma_ratio"] = (df["CV"] + eps) / (df["Magallana"] + eps)

# Filter expansions/contractions
df_focus = df[df["Category"].isin(["CV_expanded", "CV_contracted"])].copy()

# Top N for each
TOP_N = 20
top_exp = df_focus[df_focus["Category"] == "CV_expanded"].sort_values(
    "cv_vs_mgma_ratio", ascending=False
).head(TOP_N)

top_contr = df_focus[df_focus["Category"] == "CV_contracted"].sort_values(
    "cv_vs_mgma_ratio", ascending=True
).head(TOP_N)

# Label for y-axis
top_exp["Label"] = top_exp["Annotation"] + " | " + top_exp["Orthogroup"]
top_contr["Label"] = top_contr["Annotation"] + " | " + top_contr["Orthogroup"]

# Concatenate with preserved order: expanded (top) → contracted (bottom)
df_top = pd.concat([top_exp, top_contr], axis=0)

# Extract heatmap data and set Y-axis order
heatmap_df = df_top.set_index("Label")[["CV", "MG", "MA","OE"]]
heatmap_df = heatmap_df.loc[df_top["Label"].values[::-1]]  # bottom: contractions, top: expansions

print(heatmap_df)
heatmap_df.to_csv("topOG_cv_vs_maga.csv")
# Plot
plt.figure(figsize=(8, 10))
sns.heatmap(
    heatmap_df,
    annot=False,
    fmt="d",
    cmap="YlGnBu",
    linewidths=0.5,
    linecolor='gray',
    cbar_kws={"label": "Gene count","shrink":0.3}

)

# plt.title("CV vs Magallana Orthogroups: Expansion (top) and Contraction (bottom)")
plt.xlabel("Species")
plt.ylabel("")
plt.tight_layout()
plt.savefig("heatmap_cv_expansion_contraction_ordered.png", dpi=300)
plt.show()


