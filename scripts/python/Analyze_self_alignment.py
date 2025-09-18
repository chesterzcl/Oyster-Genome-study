import pandas as pd
import re
import os
from intervaltree import Interval, IntervalTree

# Set working directory if needed
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

coords_file = "cv_vs_cv_all.filtered.coords"

columns = [
    "ref_start", "ref_end",
    "qry_start", "qry_end",
    "ref_len", "qry_len",
    "identity", 
    "ref_strand", "qry_strand",
    "ref_chr", "qry_chr"
]

df = pd.read_csv(coords_file, sep="\t", names=columns)

# === Step 2: Convert numeric columns ===
numeric_cols = [
    "ref_start", "ref_end", "qry_start", "qry_end",
    "ref_len", "qry_len", "identity"
]
df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")
df.dropna(inplace=True)

# === Step 3: Remove exact self-alignments ===
df = df[~(
    (df["ref_chr"] == df["qry_chr"]) &
    (df["ref_start"] == df["qry_start"]) &
    (df["ref_end"] == df["qry_end"]) &
    (df["ref_len"] == df["qry_len"])
)]

# === Step 4: Compute alignment metrics ===
df["alignment_length"] = (df["ref_end"] - df["ref_start"]).abs()
df["ref_mid"] = (df["ref_start"] + df["ref_end"]) // 2
df["qry_mid"] = (df["qry_start"] + df["qry_end"]) // 2
df["mid_distance"] = (df["ref_mid"] - df["qry_mid"]).abs()

# === Step 4.5: Remove overlapping ref/qry intervals on the same chromosome ===
df = df[
    (df["ref_chr"] != df["qry_chr"]) |
    (df["ref_end"] < df["qry_start"]) |
    (df["qry_end"] < df["ref_start"])
]

# === Step 5: Remove overlapping alignments ===
def remove_overlapping_alignments(df):
    result_rows = []

    for (ref_chr, qry_chr), group in df.groupby(["ref_chr", "qry_chr"]):
        group = group.sort_values(by=["identity", "alignment_length"], ascending=[False, False])

        ref_tree = IntervalTree()
        qry_tree = IntervalTree()

        for _, row in group.iterrows():
            ref_interval = (min(int(row["ref_start"]), int(row["ref_end"])),
                            max(int(row["ref_start"]), int(row["ref_end"])))
            qry_interval = (min(int(row["qry_start"]), int(row["qry_end"])),
                            max(int(row["qry_start"]), int(row["qry_end"])))

            if not ref_tree.overlaps(*ref_interval) and not qry_tree.overlaps(*qry_interval):
                result_rows.append(row)
                ref_tree.add(Interval(*ref_interval))
                qry_tree.add(Interval(*qry_interval))

    return pd.DataFrame(result_rows)

df = remove_overlapping_alignments(df)

# Generate canonical keys to identify flipped duplicates
def get_canonical_key(row):
    ref = (row["ref_chr"], min(row["ref_start"], row["ref_end"]), max(row["ref_start"], row["ref_end"]))
    qry = (row["qry_chr"], min(row["qry_start"], row["qry_end"]), max(row["qry_start"], row["qry_end"]))
    return tuple(sorted([ref, qry]))

df["dup_key"] = df.apply(get_canonical_key, axis=1)

# Drop duplicate symmetric matches
df = df.drop_duplicates(subset="dup_key").drop(columns=["dup_key"])

# === Step 6: Classify duplication types ===
def classify(row):
    if row["ref_chr"] == row["qry_chr"]:
        if row["mid_distance"] <= 100_000:
            return "Tandem Duplication"
        else:
            return "Segmental Duplication"
    else:
        return "Segmental Duplication (Inter-chromosomal)"

df["category"] = df.apply(classify, axis=1)

# === Step 7: Flag large SV-like candidates ===
df["sv_candidate"] = (df["alignment_length"] >= 10_000) & (df["identity"] >= 95)

# === Step 8: Save outputs ===
df.to_csv("coords_classified.tsv", sep="\t", index=False)
df[df["category"] == "Tandem Duplication"].to_csv("coords_tandem.tsv", sep="\t", index=False)
df[df["category"].str.contains("Segmental")].to_csv("coords_segmental.tsv", sep="\t", index=False)
df[df["sv_candidate"]].to_csv("coords_large_sv.tsv", sep="\t", index=False)

# === Step 9: Print summary ===
print(df["category"].value_counts())
print(f"Large SV-like alignments: {df['sv_candidate'].sum()}")