import pandas as pd
from gprofiler import GProfiler

# Load your file (update path if needed)
dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/"
df = pd.read_csv(dir + "ret_vs_noret.best_filtered.annotated.valid.sorted.tsv", sep="\t", header=None)

# Assign column names
df.columns = [
    "retrogene_id", "retrogene_name", "parental_id", "parental_name", 
    "identity", "aln_len", "mismatch", "gapopen", 
    "qstart", "qend", "sstart", "send", 
    "evalue", "bitscore", "qlen", "slen", 
    "qcov", "scov"
]

# Filter out unknown or empty parental genes
df = df[~df["parental_id"].isin(["-", "", "Unknown_gene", "unknown_gene"])]

print(df)

# --- group by parental gene name --------------------------------------
grouped = (
    df.groupby("parental_name")                # ← the key change
      .agg(
          n_rows             = ("retrogene_id", "size"),          # total BLAST rows
          n_retrogenes       = ("retrogene_id", "nunique"),       # unique retro-IDs
          retro_ids          = ("retrogene_id",
                               lambda x: ";".join(sorted(set(x)))),  # list of retro-IDs
          mean_identity      = ("identity", "mean"),
          mean_query_cov     = ("qcov", "mean"),
          mean_subject_cov   = ("scov", "mean")
      )
      .reset_index()
      .sort_values("n_retrogenes", ascending=False)
)

print(grouped.head())          # quick sanity-check

grouped.to_csv(dir+ "retrogene_parental_groups.tsv", sep="\t", index=False)
# Extract unique list of parental gene symbols (assumes these are human orthologs)
parent_genes = df["parental_name"].dropna().unique().tolist()

# Run GO enrichment
gp = GProfiler(return_dataframe=True)
go_results = gp.profile(
    organism="hsapiens",  # Adjust if using another species
    query=parent_genes,
    no_evidences=False
)

# Filter to only GO:BP terms and p < 0.05
go_bp = go_results[(go_results["source"] == "GO:BP") & (go_results["p_value"] < 0.05)]

# Sort by p-value and save
go_bp = go_bp.sort_values("p_value")
go_bp.to_csv(dir + "retrogene_parent_GO_BP_enrichment.tsv", sep="\t", index=False)

print(list(go_bp))

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Sort by p-value and select top N terms
top_n = 10
top_go = go_bp.sort_values("p_value").head(top_n)

# Plot
plt.figure(figsize=(10, 6))
bars = plt.barh(
    y=top_go["name"][::-1], 
    width=-top_go["p_value"].apply(lambda p: np.log10(p))[::-1], 
    color="steelblue"
)

plt.xlabel("-log10(p-value)")
plt.title(f"Top {top_n} Enriched GO:Biological Process Terms")
plt.tight_layout()
plt.grid(axis="x", linestyle="--", alpha=0.5)
plt.show()
