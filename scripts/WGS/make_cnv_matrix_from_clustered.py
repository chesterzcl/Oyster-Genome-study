import pandas as pd

# === Load BED ===
df = pd.read_csv("clustered_beds/all_samples_clustered.bed", sep="\t", header=None,
                 names=["chr", "start", "end", "type", "sample", "cluster"])

# === Create unique CNV cluster ID ===
df["cluster_id"] = (df["chr"].astype(str) + ":" +
                    df["start"].astype(str) + "-" +
                    df["end"].astype(str) + "_" +
                    df["type"] + "_cl" + df["cluster"].astype(str))

# === Build presence matrix ===
presence = df.pivot_table(index="cluster_id", columns="sample", aggfunc="size", fill_value=0)

# === Optional: sort by chromosome and start ===
presence["sort_chr"] = presence.index.str.extract(r'^(.*?):')[0]
presence["sort_start"] = presence.index.str.extract(r':(\d+)-')[0].astype(int)
presence = presence.sort_values(by=["sort_chr", "sort_start"]).drop(columns=["sort_chr", "sort_start"])

# === Output ===
presence.to_csv("cnv_cluster_presence_matrix.tsv", sep="\t")
presence.to_csv("cnv_cluster_presence_matrix.csv")

# === Optional: per-sample CNV counts ===
presence.sum(axis=0).to_csv("cnv_counts_per_sample.tsv", sep="\t")

summary = presence.sum(axis=0)
summary.to_csv("cnv_counts_per_sample.tsv", sep="\t")

