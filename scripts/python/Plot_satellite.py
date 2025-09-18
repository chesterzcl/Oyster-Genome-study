import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load filtered BED file
bed_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/satellite_analysis_satellite_filtered.bed"
fai_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/primary_dedup_chr.fa.fai"
out_dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/"

window_size = 100000

# Load satellite BED
df = pd.read_csv(bed_file, delim_whitespace=True, header=None, names=[
    "chrom", "start", "end", "motif", "score", "strand",
    "period", "copies", "match", "indel", "length"
])



# Load genome index (.fai)
fai = pd.read_csv(fai_file,sep="\t",header=None,usecols=[0,1],names=["chrom", "length"])

# Basic stats
total_arrays = len(df)
total_bp = df["length"].sum()
mean_len = df["length"].mean()
median_len = df["length"].median()
most_common_motif = df["motif"].value_counts().idxmax()
motif_counts = df["motif"].value_counts()

# Period size distribution
period_stats = df["period"].describe()

# Scaffold summary
scaffold_bp = df.groupby("chrom")["length"].sum().sort_values(ascending=False)

# Print summary
print(f"Total arrays: {total_arrays}")
print(f"Total satellite bp: {total_bp}")
print(f"Mean array length: {mean_len:.1f}")
print(f"Median array length: {median_len}")
print(f"Most common motif: {most_common_motif}")
print("\nTop 5 motifs:\n", motif_counts.head(40))
print("\nPeriod size distribution:\n", period_stats)
print("\nTop 5 scaffolds by satellite content:\n",scaffold_bp.head(10))

chrom_summary = df.groupby("chrom").agg(
    satellite_count=("length", "count"),
    total_satellite_bp=("length", "sum"),
    mean_satellite_length=("length", "mean"),
    median_satellite_length=("length", "median")
).sort_values(by="total_satellite_bp", ascending=False)

print("Satellite Distribution by Chromosome:")
chrom_summary = chrom_summary.sort_values(by="chrom")
chrom_summary = chrom_summary.merge(fai, on="chrom")
chrom_summary["satellite_fraction"]=chrom_summary["total_satellite_bp"]/chrom_summary["length"]*100
print(chrom_summary.head(10))
chrom_summary.to_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/satellite_by_chromosome.txt", sep="\t", float_format="%.3f")

# Bin genome and count satellite coverage
stlt_seq1="TCACCCAAGGATGCTTTGTGCCAAGTTTGGTTGAAATTGGCTCAGTGGTTCTGGAGAAGAAGATTTTTAAATTTCGTCAATGTATTTTCACCATTTCGTAATTATCTCCCCTTGGAAAAGGGCGGGGTCCTTCATTTGAACAAACTTGAATCCCCT"
stlt_seq2="ATCTTCTTCTCCAGAACCACTGAGCCAATTTCAACCAAACTTGGCACAAAGCATCCTTGGGTGAAGGGGATTCAAGTTTGTTCAAATGAAGGACCCCGCCCTTTTCCAAGGGGAGATAATTACGAAATAGTGAAAATACATTGACGAAATTTAAAA"

tel_seq1="CCCTAA"
tel_seq2="TTAGGG"

def get_rotations(sequence: str) -> set:
    return {sequence[i:] + sequence[:i] for i in range(len(sequence))}

tel_set1=get_rotations(tel_seq1)
tel_set2=get_rotations(tel_seq2)

stlt_set1=get_rotations(stlt_seq1)
stlt_set2=get_rotations(stlt_seq2)

# print(df[df["motif"].str.contains(seq1, regex=True)])
# print(df[df["motif"].str.contains(seq2, regex=True)])

tel_set_combined=tel_set1|tel_set2
stlt_set_combined=stlt_set1|stlt_set2

# Filter once for all reverse-complement motif instances
df_filtered = df[df["motif"].isin(tel_set_combined)].copy()
print(df_filtered)

# Plot satellite density across chromosomes (no motif filtering)
for chrom, length in fai.values:
    bins = np.arange(0, length + window_size, window_size)
    bin_labels = bins[:-1]

    sat = df[df["chrom"] == chrom].copy()
    if sat.empty:
        continue

    chr_id = chrom.split('_')[-1]
    sat["bin"] = pd.cut(sat["start"], bins=bins, labels=bin_labels, right=False)
    binned = sat.groupby("bin")["length"].sum().reindex(bin_labels, fill_value=0).reset_index()
    binned.columns = ["start", "sat_bp"]
    binned["start"] = binned["start"].astype(int)

    # Plot with fixed y-axis
    plt.figure(figsize=(10, 3))
    plt.plot(binned["start"] / 1e6, binned["sat_bp"], lw=1, color="steelblue")
    plt.title(f"Chr{chr_id} – Total satellite density")
    plt.xlabel("Position (Mb)")
    plt.ylabel("Satellite bp")
    plt.ylim(0, 200000)  # Adjust this based on your genome size and density
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Chr{chr_id}_total_satellite_density.png", dpi=300)
    plt.close()


# Plot by chromosome
for chrom, length in fai.values:
    bins = np.arange(0, length + window_size, window_size)
    bin_labels = bins[:-1]

    sat = df_filtered[df_filtered["chrom"] == chrom].copy()
    if sat.empty:
        continue

    chr_id = chrom.split('_')[-1]
    sat["bin"] = pd.cut(sat["start"], bins=bins, labels=bin_labels, right=False)
    binned = sat.groupby("bin")["length"].sum().reindex(bin_labels, fill_value=0).reset_index()
    binned.columns = ["start", "sat_bp"]
    binned["start"] = binned["start"].astype(int)

    # Plot single-color line
    plt.figure(figsize=(10, 3))
    plt.plot(binned["start"] / 1e6, binned["sat_bp"], lw=1, color="darkorange")
    plt.title(f"Chr{chr_id} – Telomeric satellite density")
    plt.xlabel("Position (Mb)")
    plt.ylabel("Satellite bp")
    plt.ylim(0, 8000)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Chr{chr_id}_telo_density.png", dpi=300)
    plt.close()

df_filtered.to_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/telo_stlt.txt", sep="\t", float_format="%.3f")

df_filtered = df[df["motif"].isin(stlt_set_combined)].copy()
df_filtered.to_csv("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/nuc_stlt.txt", sep="\t", float_format="%.3f")
print(df_filtered)
# Plot by chromosome
for chrom, length in fai.values:
    bins = np.arange(0, length + window_size, window_size)
    bin_labels = bins[:-1]

    sat = df_filtered[df_filtered["chrom"] == chrom].copy()
    if sat.empty:
        continue

    chr_id = chrom.split('_')[-1]
    sat["bin"] = pd.cut(sat["start"], bins=bins, labels=bin_labels, right=False)
    binned = sat.groupby("bin")["length"].sum().reindex(bin_labels, fill_value=0).reset_index()
    binned.columns = ["start", "sat_bp"]
    binned["start"] = binned["start"].astype(int)

    # Plot single-color line
    plt.figure(figsize=(10, 3))
    plt.plot(binned["start"] / 1e6, binned["sat_bp"], lw=1, color="darkorange")
    plt.title(f"Chr{chr_id} – Nucleosome wrapping satellite density")
    plt.xlabel("Position (Mb)")
    plt.ylabel("Satellite bp")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Chr{chr_id}_stlt_density.png", dpi=300)
    plt.close()

