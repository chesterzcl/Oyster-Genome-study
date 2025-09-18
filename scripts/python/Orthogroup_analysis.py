import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore

# os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/Orthogroups")
# os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final/Orthogroups_CV_MG")
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/Orthogroups/")

# === Step 1: Parse Orthogroups.tsv to map Orthogroup → CV gene IDs ===
def parse_orthogroups(file):
    ortho2cv = {}
    with open(file) as f:
        for line in f:
            if line.startswith("Orthogroup"): continue
            parts = line.strip().split("\t")
            ortho_id = parts[0]
            if len(parts) > 1:
                cv_field = parts[1]
                cv_genes = [x.split("|")[-1] for x in cv_field.split(", ") if x.startswith("Cvir|")]
                ortho2cv[ortho_id] = cv_genes
    return ortho2cv

# === Step 2: Parse GFF for gene ID → gene name ===
def parse_gff_for_gene_names(gff_file):
    gene_to_name = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            feature, attr_str = parts[2], parts[8]
            if feature not in ["gene", "mRNA"]: continue
            attrs = dict(field.split("=", 1) for field in attr_str.split(";") if "=" in field)
            gene_id = attrs.get("ID", "").split(".")[0]
            gene_name = attrs.get("Name")
            if gene_id and gene_name:
                gene_to_name[gene_id] = gene_name
    return gene_to_name

# === Step 3: Heatmap plotting ===
def plot_heatmap(df_sub, species_cols, filename, title, ortho2cv, gene_name_map, figsize=(10, 12), cmap="vlag"):
    mat = df_sub[species_cols].apply(pd.to_numeric, errors="coerce").pipe(lambda x: np.log2(x + 1))
    z_mat = mat.apply(lambda row: zscore(row.values, nan_policy="omit"), axis=1, result_type="expand")
    z_mat.columns = species_cols

    # Generate improved y-axis labels
    labels = []
    for ortho in df_sub["Orthogroup"]:
        gene_ids = ortho2cv.get(ortho, [])
        names = [gene_name_map.get(gid) for gid in gene_ids if gene_name_map.get(gid)]
        label = names[0] if names else ortho
        labels.append(label)

    z_mat.index = labels

    plt.figure(figsize=figsize)
    sns.heatmap(
        z_mat,
        cmap=cmap,
        linewidths=0.3,
        linecolor="#E0E0E0",
        cbar_kws={"label": "Row z-score (log₂ copy-number)"},
        xticklabels=True,
        yticklabels=True
    )
    plt.title(title, fontsize=14, pad=12)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()

# === Step 4: Load all data ===
ortho2cv = parse_orthogroups("Orthogroups.tsv")
gene_name_map = parse_gff_for_gene_names("CV.gff")

df = pd.read_csv("Orthogroups.GeneCount.annotated.tsv", sep="\t")
df = df.rename(columns={
    "C.virginica": "CV",
    "M.angulata" : "MA",
    "M.edulis"   : "ME",
    "M.gigas"    : "MG",
    "O.edulis":"OE"
})

print(df)
# === Step 5: Fold-change and filtering ===
df["mean_MA_MG"] = (df["MA"] + df["MG"]) / 2
df["mean_oyster"] = (df["CV"] + df["MA"] + df["MG"]) / 3
df["fold_CV_vs_MA_MG"] = df["CV"] / (df["mean_MA_MG"] + 1e-5)
df["fold_oysters_vs_ME"] = df["mean_oyster"] / (df["ME"] + 1e-5)

fc_cutoff = 2

cv_enriched = df[
    (df[["CV", "MA", "MG"]] >= 1).all(axis=1) &
    ((df["fold_CV_vs_MA_MG"] >= fc_cutoff) |
     (df["fold_CV_vs_MA_MG"] <= 1 / fc_cutoff))
]

oyster_expanded = df[
    (df[["CV", "MA", "MG", "ME"]] > 0).all(axis=1) &
    ((df["fold_oysters_vs_ME"] >= fc_cutoff) |
     (df["fold_oysters_vs_ME"] <= 1 / fc_cutoff))
]

# === Step 6: Filter only orthogroups with named CV genes ===
def has_named_cv_gene(ortho_id):
    gene_ids = ortho2cv.get(ortho_id, [])
    return any(gid in gene_name_map for gid in gene_ids)

cv_enriched = cv_enriched[cv_enriched["Orthogroup"].apply(has_named_cv_gene)]
oyster_expanded = oyster_expanded[oyster_expanded["Orthogroup"].apply(has_named_cv_gene)]

cv_extremes = pd.concat([
    cv_enriched.nlargest(20, "fold_CV_vs_MA_MG"),
    cv_enriched.nsmallest(20, "fold_CV_vs_MA_MG")
])

oy_extremes = pd.concat([
    oyster_expanded.nlargest(20, "fold_oysters_vs_ME"),
    oyster_expanded.nsmallest(20, "fold_oysters_vs_ME")
])

# === Step 7: Final plots ===
plot_heatmap(
    cv_extremes,
    ["CV", "MA", "MG"],
    "heatmap_CV_enriched_labeled.png",
    "CV vs MA+MG – extreme orthogroups",
    ortho2cv, gene_name_map,
    figsize=(10, 12),
    cmap="vlag"
)

plot_heatmap(
    oy_extremes,
    ["CV", "MA", "MG", "ME"],
    "heatmap_Oyster_vs_ME_labeled.png",
    "Oysters vs M. edulis – extreme orthogroups",
    ortho2cv, gene_name_map,
    figsize=(11, 13),
    cmap="vlag"
)

print("Heatmaps saved with only named CV gene orthogroups.")



# Get orthogroup IDs with at least one named CV gene
valid_orthogroups = []
for ortho, genes in ortho2cv.items():
    for g in genes:
        if g in gene_name_map:
            valid_orthogroups.append(ortho)
            break

# Filter cv_extremes and oy_extremes based on this
cv_extremes = cv_extremes[cv_extremes["Orthogroup"].isin(valid_orthogroups)]
oy_extremes = oy_extremes[oy_extremes["Orthogroup"].isin(valid_orthogroups)]

# For CV-enriched set
cv_genes_named = []
for ortho in cv_extremes["Orthogroup"]:
    for gid in ortho2cv.get(ortho, []):
        if gid in gene_name_map:
            cv_genes_named.append(gene_name_map[gid])

# For oyster-expanded set
oy_genes_named = []
for ortho in oy_extremes["Orthogroup"]:
    for gid in ortho2cv.get(ortho, []):
        if gid in gene_name_map:
            oy_genes_named.append(gene_name_map[gid])

from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

def run_go_and_plot(gene_list, list_name, top_n=20,
                    organism="hsapiens", out_dir=".", palette="Greens_d"):
    """
    gene_list  – list of gene symbols
    list_name  – basename for output files (prefix)
    """
    if not gene_list:
        print(f"[WARN] {list_name}: empty gene list – skipping.")
        return

    # --- 1. Run enrichment -------------------------------
    res = gp.profile(
        organism=organism,
        query=gene_list,
        sources=["GO:BP"],
        no_evidences=False
    )

    if res.empty:
        print(f"[WARN] {list_name}: no GO terms returned.")
        return

    # --- 2. Save full GO results --------------------------
    tsv_path = os.path.join(out_dir, f"{list_name}_GO_BP.tsv")
    res.to_csv(tsv_path, sep="\t", index=False)
    print(f"GO table saved  →  {tsv_path}")

    # --- 3. Filter for Biological Process -----------------
    bp = res[res["source"] == "GO:BP"]
    if bp.empty:
        print(f"[WARN] {list_name}: no GO:BP terms found.")
        return

    # --- 4. Select top_n by smallest p-value --------------
    top = bp.nsmallest(top_n, "p_value").copy()
    top["-log10(padj)"] = -np.log10(top["p_value"])

    # --- 5. Plot -------------------------------------------
    plt.figure(figsize=(8, 0.5 * len(top) + 1))  # Adjust height to fit labels
    sns.barplot(
        data=top,
        x="-log10(padj)",
        y="name",
        hue="source",
        palette=palette
    )
    plt.xlabel("-log10(FDR)", fontsize=11)
    plt.ylabel("GO Biological Process", fontsize=11)
    plt.title(f"Top {top_n} GO:BP terms – {list_name}", fontsize=12)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=8)
    plt.legend(loc="lower right", fontsize=8)
    plt.tight_layout()

    png_path = os.path.join(out_dir, f"{list_name}_GO_BP.png")
    plt.savefig(png_path, dpi=300)
    plt.close()
    print(f"GO plot saved → {png_path}")

# -------------------------------------------------------------------------
#  Run for the two filtered gene sets you already built
# -------------------------------------------------------------------------
run_go_and_plot(cv_genes_named, "CV_enriched_BP",  top_n=30, palette="Blues_d")
run_go_and_plot(oy_genes_named, "Oyster_vs_ME_BP", top_n=30, palette="Oranges_d")


from gseapy import prerank

def run_gsea_analysis(df, ortho2cv, gene_name_map, fc_column, output_prefix, organism="Human", outdir="gsea_results", min_size=3, max_size=1000):
    # Step 1: Construct ranked list
    ranked = []
    for _, row in df.iterrows():
        ortho = row["Orthogroup"]
        score = row[fc_column]
        gene_ids = ortho2cv.get(ortho, [])
        for gid in gene_ids:
            if gid in gene_name_map:
                ranked.append((gene_name_map[gid], score))
    
    ranked_df = pd.DataFrame(ranked, columns=["gene", "score"])
    ranked_df = ranked_df.groupby("gene").mean().sort_values("score", ascending=False)

    # Output file
    rnk_file = f"{output_prefix}.rnk"
    ranked_df.to_csv(rnk_file, sep="\t", header=False)
    print(f"Ranked gene list saved → {rnk_file}")

    # Step 2: Run GSEA pre-ranked
    prerank(
        rnk=rnk_file,
        gene_sets="GO_Biological_Process_2021",
        organism=organism,
        outdir=os.path.join(outdir, output_prefix),
        format="png",
        min_size=min_size,
        max_size=max_size,
        permutation_num=100,  # Reduce for faster testing
        seed=42,
        no_plot=False,
        verbose=True
    )

    print(f"GSEA completed → {os.path.join(outdir, output_prefix)}")


# run_gsea_analysis(df, ortho2cv, gene_name_map, "fold_CV_vs_MA_MG", "GSEA_CV_vs_MA_MG")
# run_gsea_analysis(df, ortho2cv, gene_name_map, "fold_oysters_vs_ME", "GSEA_Oysters_vs_ME")

