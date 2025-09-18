import pandas as pd

# === CONFIGURATION ===
bed_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/repeat_analysis/all_repeats.bed"
out_dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/repeat_analysis/"
genome_size = 532076140

# === LOAD BED FILE ===
df = pd.read_csv(bed_file, sep="\t", names=[
    "chrom", "start", "end", "repeat_id", "dot", "class/family"
])

# === PARSE CLASSES & LENGTH ===
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)
df["length"] = df["end"] - df["start"]
df["class"] = df["class/family"].str.extract(r"^([^/]+)")
df["subclass"] = df["class/family"].str.extract(r"/(.+)$")

# === SUMMARY 1: BY CLASS ===
class_summary = (
    df.groupby("class")["length"]
    .sum()
    .reset_index()
    .rename(columns={"length": "total_bp"})
)
class_summary["percent_genome"] = (class_summary["total_bp"] / genome_size) * 100
class_summary = class_summary.sort_values(by="total_bp", ascending=False)

# === SUMMARY 2: BY CLASS + SUBCLASS ===
subclass_summary = (
    df.groupby(["class", "subclass"])["length"]
    .sum()
    .reset_index()
    .rename(columns={"length": "total_bp"})
)
subclass_summary["percent_genome"] = subclass_summary["total_bp"] / genome_size * 100
subclass_summary = subclass_summary.sort_values(by=["class", "total_bp"], ascending=[True, False])

# === OPTIONAL: ADD RANK PER CLASS ===
subclass_summary["rank_within_class"] = (
    subclass_summary.groupby("class")["total_bp"]
    .rank(method="dense", ascending=False).astype(int)
)

# === EXPORT BOTH TO TXT ===
class_summary.to_csv(out_dir+"repeat_class_summary.txt", sep="\t", index=False)
subclass_summary.to_csv(out_dir+"repeat_subclass_summary.txt", sep="\t", index=False)