import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


genes=pd.read_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/gene_analysis/genes.bed",sep="\t",header=None,names=["chrom","start","end","gene_id"])
eggnog=pd.read_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/eggnog_annotation.emapper.annotations",sep="\t",comment="#")


total_queries = eggnog.shape[0]
with_function = eggnog["Description"].notna().sum()

print(f"Total predicted proteins: {total_queries}")
print(f"With functional description: {with_function} ({with_function / total_queries:.1%})")


top_genes = eggnog["Preferred_name"].value_counts().head(10)
top_desc = eggnog["Description"].value_counts().head(10)


print("Top Preferred Names:\n", top_genes)
print("\nTop Descriptions:\n", top_desc)


has_name = eggnog[
    eggnog["Preferred_name"].notna() &
    (eggnog["Preferred_name"].str.strip() != "-")
].copy()


has_name["gene_id"] = has_name["query"].str.extract(r"^(g\d+)")
unique_genes = has_name["gene_id"].unique()

print(f"Found {len(has_name)} annotated transcripts, representing {len(unique_genes)} genes.")
print("Sample entries:")
print(has_name[["query","gene_id","Preferred_name","Description"]].head())


cog_string = eggnog["COG_category"].dropna().str.cat().replace("-", "")
cog_counts = pd.Series(list(cog_string)).value_counts()

cog_labels = {
    "J": "Translation, ribosomal structure",
    "K": "Transcription",
    "L": "Replication, recombination",
    "M": "Cell wall/membrane",
    "N": "Cell motility",
    "O": "Post-translational modification",
    "T": "Signal transduction",
    "U": "Intracellular trafficking",
    "V": "Defense mechanisms",
    "X": "Mobilome (transposons, prophages)",
    "R": "General function",
    "S": "Unknown function",
    # ... (you can add more if needed)
}

# Annotate categories
cog_df=pd.DataFrame({
    "Category": cog_counts.index,
    "Count": cog_counts.values,
    "Function": [cog_labels.get(c, "Other") for c in cog_counts.index]
})

print(cog_df.head(10))


te_keywords = ["transposase", "retrotransposon", "reverse transcriptase", "gag", "pol", "integrase", "mobile element"]
te_related = eggnog["Description"].str.lower().str.contains("|".join(te_keywords), na=False)

print(f"TE-related annotations: {te_related.sum()} ({te_related.mean():.1%})")


# Rename the relevant ID column
eggnog=eggnog.rename(columns={"query":"transcript_id"})


# Extract gene ID (e.g., g1.t1 → g1)
eggnog["gene_id"]=eggnog["transcript_id"].str.extract(r"^(g\d+)")

# === Define high-confidence protein coding genes ===
protein_coding = eggnog[
    (eggnog["Description"].notna()) |
    (eggnog["COG_category"].notna()) |
    (eggnog["KEGG_ko"].notna()) |
    (eggnog["PFAMs"].notna())
].copy()

annotated = pd.merge(genes,protein_coding[[
    "gene_id", "Preferred_name", "Description", "COG_category", "KEGG_ko","PFAMs"
]], on="gene_id", how="left")

# === Define annotation type (protein_coding / unknown) ===
annotated["feature_type"] = annotated["Description"].notna().map({True:"protein_coding",False: "unknown_gene"})


# Or for genome browser (name = function)
annotated[["chrom", "start", "end", "gene_id", "Preferred_name", "feature_type"]].drop_duplicates().reset_index(drop=True).to_csv(
    "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/genes_function.bed", sep="\t", header=False, index=False)
print(annotated["feature_type"])



