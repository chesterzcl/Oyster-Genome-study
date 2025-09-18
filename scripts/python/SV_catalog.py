import pysam
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === Path to your VCF file ===
dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/sv/"
vcf_path = dir+"merged_sv.sorted.vcf"

# === Open the VCF ===
vcf = pysam.VariantFile(vcf_path)

# === Collect data ===
records = []

for rec in vcf.fetch():
    info = rec.info
    svtype = info.get("SVTYPE", "NA")
    svlen = abs(info.get("SVLEN", 0))  # abs in case of negative DEL
    supp = int(info.get("SUPP", 0))
    qual = rec.qual if rec.qual else 0.0

    records.append({
        "CHROM": rec.chrom,
        "POS": rec.pos,
        "SVTYPE": svtype,
        "SVLEN": svlen,
        "SUPP": supp,
        "QUAL": qual
    })

# === Convert to DataFrame ===
df = pd.DataFrame(records)

# === Summary statistics ===
print("\n=== Count by SVTYPE ===")
print(df["SVTYPE"].value_counts())

print("\n=== Summary by SVLEN ===")
print(df["SVLEN"].describe())

print("\n=== Summary by QUAL ===")
print(df["QUAL"].describe())

print("\n=== Summary by SUPP (supporting samples) ===")
print(df["SUPP"].describe())

# === Optional: Save as CSV ===
df.to_csv("sv_summary.csv", index=False)
print("\nSaved summary to sv_summary.csv")

# Filter reliable SVs
filtered_df = df[
    (df["SVLEN"].abs() >= 50) &
    (df["SVLEN"].abs() <= 1e5) &
    (df["QUAL"] >= 500) &
    (df["SUPP"] >= 2) &
    (df["SUPP"] <= 18) &
    (df["SVTYPE"].isin(["DEL", "DUP", "INS", "INV"]))
].copy()

print(f"Filtered SVs retained: {len(filtered_df)} / {len(df)}")

# Optionally export
filtered_df.to_csv(dir+"mantle_sv_filtered.csv", index=False)

df=filtered_df

# Set Seaborn style
sns.set(style="whitegrid")

# Set style globally
sns.set_context("paper", font_scale=1.6)
sns.set_style("whitegrid")


# Add a new column for log10(SVLEN)
df["LOG10_SVLEN"] = np.log10(df["SVLEN"].replace(0, np.nan)).dropna()


# === 1. SV Type Frequency ===
plt.figure(figsize=(7, 5))
ax = sns.countplot(data=df, x="SVTYPE", order=df["SVTYPE"].value_counts().index, palette="muted")
ax.set_title("SV Type Distribution")
ax.set_xlabel("SV Type")
ax.set_ylabel("Count")
sns.despine()
plt.tight_layout()
plt.savefig(dir + "sv_type_frequency.png")  # Use vector format
plt.close()

# === Log10(SVLEN) Histogram ===
plt.figure(figsize=(7, 5))
sns.histplot(data=df, x="LOG10_SVLEN", bins=50, kde=True, color="steelblue")
plt.title("SV Length Distribution")
plt.xlabel("log10(SV Length)")
plt.ylabel("Count")
sns.despine()
plt.tight_layout()
plt.savefig(dir + "sv_length_log10_distribution.png")
plt.close()

# === 3. QUAL Distribution ===
plt.figure(figsize=(7, 5))
ax = sns.histplot(data=df, x="QUAL", bins=50, kde=True, color="darkgreen")
ax.set_title("SV QUAL Score Distribution")
ax.set_xlabel("QUAL Score")
ax.set_ylabel("Count")
sns.despine()
plt.tight_layout()
plt.savefig(dir + "sv_qual_distribution.png")
plt.close()

# === 4. Boxplot: SV Length by SV Type (<100kb) ===
plt.figure(figsize=(8, 5))
filtered = df[df["SVLEN"] < 100000]
filtered["LOG10_SVLEN"] = np.log10(filtered["SVLEN"].replace(0, np.nan))

# === Boxplot of log10 SVLEN by SVTYPE ===
flier_props = dict(marker='.', markersize=5, linestyle='none', markeredgecolor='gray')

plt.figure(figsize=(8, 5))
sns.boxplot(
    data=df,
    x="SVTYPE",
    y="LOG10_SVLEN",
    palette="Set2",
    flierprops=flier_props
)
plt.title("SV Length by SV Type")
plt.xlabel("SV Type")
plt.ylabel("log10(SV Length)")
sns.despine()
plt.tight_layout()
plt.savefig(dir + "sv_length_log10_by_type.png")
plt.close()

print("Saved all publication-grade SV distribution plots.")


# === Rewind and re-open original VCF to re-iterate over full records ===
vcf = pysam.VariantFile(vcf_path)

# === Prepare VCF writer ===
filtered_vcf_path = dir + "mantle_sv_filtered.vcf"
filtered_vcf = pysam.VariantFile(filtered_vcf_path, 'w', header=vcf.header)

# === Reapply filter and write passing records ===
count_written = 0
for rec in vcf.fetch():
    info = rec.info
    svtype = info.get("SVTYPE", "NA")
    svlen = abs(info.get("SVLEN", 0))
    supp = int(info.get("SUPP", 0))
    qual = rec.qual if rec.qual else 0.0

    if (
        50 <= svlen <= 1e5 and
        qual >= 500 and
        2<= supp <=18 and
        svtype in ["DEL", "DUP", "INS", "INV"]
    ):
        filtered_vcf.write(rec)
        count_written += 1

filtered_vcf.close()
print(f"Saved filtered VCF with {count_written} variants to: {filtered_vcf_path}")