import pandas as pd
import numpy as np
from collections import defaultdict
import os
import pyranges as pr


os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")
# === Parameters ===
repeat_bed = "all_repeats.bed"         # Input repeat BED with class in 4th column
fai_file = "primary_dedup_chr_masked_hp_sealed.fa.fai"             # Input .fai file with chromosome sizes
window_size = 200000                     # Bin size in bp
output_file = f"repeat_TE_stacked_{window_size}.tsv" # Output file

# 1. Load repeats and keep four subclasses
focus_classes = ["DNA", "LTR", "LINE", "SINE"]
# focus_classes = ["DNA", "LTR", "LINE", "SINE","Unknown"]
df = pd.read_csv(repeat_bed, sep="\t", header=None,
                 names=["chrom", "start", "end", "id", "dot","strand", "class/family"])
df["class"] = df["class/family"].str.extract(r"^([^/]+)")     # take text before '/'
df = df[df["class"].isin(focus_classes)]
df["length"] = df["end"].astype(int) - df["start"].astype(int)
df["chrom"]  = df["chrom"].astype(str).str.strip()

# 2. Load chromosome sizes and build bin table
chrom_sizes = pd.read_csv(fai_file, sep="\t", header=None,
                          usecols=[0, 1], names=["chrom", "size"])
chrom_sizes["chrom"] = chrom_sizes["chrom"].astype(str).str.strip()

bins = []
for _, row in chrom_sizes.iterrows():
    chrom, size = row["chrom"], row["size"]
    edges = np.arange(0, size + window_size, window_size)
    bins += [
        {"chrom": chrom,
         "bin_start": int(edges[i]),
         "bin_end":   int(edges[i+1])}
        for i in range(len(edges) - 1)
    ]
bin_df = pd.DataFrame(bins)


######==========================================================================
chrom_sizes_dict = dict(zip(chrom_sizes["chrom"], chrom_sizes["size"]))

# === Compute non-overlapping repeat class lengths per chromosome ===
records = []

for (chrom, rep_class), group in df.groupby(["chrom", "class"]):
    gr = pr.PyRanges(group.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
    merged = gr.merge()
    rep_len = merged.lengths().sum()
    chrom_len = chrom_sizes_dict.get(chrom, None)
    if chrom_len:
        pct = 100 * rep_len / chrom_len
        records.append({
            "chrom": chrom,
            "repeat_class": rep_class,
            "nonoverlap_bp": rep_len,
            "chrom_length": chrom_len,
            "pct_coverage": pct
        })

# === Assemble per-chromosome percent coverage table ===
repeat_pct_df = pd.DataFrame(records)

# Pivot: rows = chromosome, columns = repeat_class, values = percentage
percent_pivot = repeat_pct_df.pivot(index="chrom", columns="repeat_class", values="pct_coverage").fillna(0)

# === Compute genome-wide repeat coverage ===
genome_records = []

for rep_class, group in df.groupby("class"):
    gr = pr.PyRanges(group.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
    merged = gr.merge()
    total_rep_len = merged.lengths().sum()
    genome_records.append({"repeat_class": rep_class, "genome_nonoverlap_bp": total_rep_len})

# Genome-wide percentage
genome_total_len = chrom_sizes["size"].sum()
genome_df = pd.DataFrame(genome_records)
genome_df["pct_coverage"] = 100 * genome_df["genome_nonoverlap_bp"] / genome_total_len
genome_row = genome_df.set_index("repeat_class")["pct_coverage"].T.to_frame().T
genome_row.index = ["Genome"]

# === Final combined table ===
final_df = pd.concat([percent_pivot, genome_row], axis=0).round(3)

# === Save or display ===
final_df.to_csv("repeat_percent_coverage.tsv", sep="\t")
######==========================================================================


# 3. Assign repeats to overlapping bins
def overlapping_bins(rep_row):
    sel = bin_df.loc[
        (bin_df["chrom"] == rep_row["chrom"]) &
        (rep_row["start"] < bin_df["bin_end"]) &
        (rep_row["end"]   > bin_df["bin_start"])
    ]
    return sel.index.tolist()

rep_to_bins = df.apply(overlapping_bins, axis=1)
mask_valid  = rep_to_bins.str.len() > 0
df = df[mask_valid].reset_index(drop=True)
rep_to_bins = rep_to_bins[mask_valid].reset_index(drop=True)

# explode repeat-bin relationships
df_exp = df.loc[df.index.repeat(rep_to_bins.str.len())].copy()
df_exp["bin_index"] = [idx for sub in rep_to_bins for idx in sub]
df_exp = df_exp.merge(bin_df, left_on="bin_index", right_index=True, suffixes=("", "_bin"))
# 4. Collapse subclasses → two groups
retro_classes = {"LINE", "SINE", "LTR"}
df_exp["group"] = np.where(df_exp["class"].isin(retro_classes),
                           "Retrotransposon", "DNA_transposon")
print(df_exp)
# 5. Sum bp per bin per group
print(df_exp.columns)
stats = (df_exp
         .groupby(["chrom", "bin_start", "bin_end", "group"])["length"]
         .sum()
         .reset_index())

# pivot to wide format
wide = (stats
        .pivot_table(index=["chrom", "bin_start", "bin_end"],
                     columns="group", values="length", fill_value=0)
        .reset_index())
wide.columns.name = None   # drop pandas pivot axis name

# ensure both columns exist
for col in ["Retrotransposon", "DNA_transposon"]:
    if col not in wide:
        wide[col] = 0

# 6. Add % coverage columns
wide["retro_pct"] = wide["Retrotransposon"] / window_size * 100
wide["dna_pct"]   = wide["DNA_transposon"] / window_size * 100

# reorder & save
out_cols = ["chrom", "bin_start", "bin_end",
            "Retrotransposon", "DNA_transposon", "retro_pct", "dna_pct"]
wide[out_cols].to_csv(output_file, sep="\t", index=False)

