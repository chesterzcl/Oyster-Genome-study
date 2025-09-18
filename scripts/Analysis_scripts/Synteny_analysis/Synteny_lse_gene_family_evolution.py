import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Optional: Set your working directory
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final")

# --- Load both files ---
df_whole = pd.read_csv("cds_unique.fa.tsv.ks.tsv", sep="\t")
df_lse = pd.read_csv("lse_900genes.fa.tsv.ks.tsv", sep="\t")

# --- Clean and standardize ---
def clean_df(df):
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["dN", "dS", "dN/dS", "node_averaged_dS_outlierexcluded"])
    df = df[(df["dN"] > 0) & (df["dS"] > 0) & (df["dN/dS"] > 0) & (df["dN"] <=10 ) & (df["dS"] <=10) & (df["dN/dS"] <=10 ) ]
    # df = df[df["node_averaged_dS_outlierexcluded"] <= 5]  # mimic WGD cap
    df["effective_aligned_bp"] = df["alignmentlength"] * df["alignmentcoverage"]
    df=df[df["effective_aligned_bp"] >= 150]
    df=df[df["alignmentidentity"] >= 0.30]
    df["Ka"] = df["dN"]
    df["Ks"] = df["node_averaged_dS_outlierexcluded"]
    # df["Ks"] = df["dS"]
    df["omega"] = df["dN/dS"]
    df["log10_Ka"] = np.log10(df["Ka"])
    df["log10_Ks"] = np.log10(df["Ks"])
    df["log10_omega"] = np.log10(df["omega"])
    return df

df_whole = clean_df(df_whole)
df_lse = clean_df(df_lse)

# --- Plotting: log10(Ks) overlay ---
plt.figure(figsize=(6, 5))
bins = np.linspace(-2.99, 1, 30)
plt.hist(df_lse["log10_Ks"], bins=bins, alpha=0.6, label="LSE", density=True, color="seagreen")
plt.hist(df_whole["log10_Ks"], bins=bins, alpha=0.6, label="Whole Genome", density=True,color="gray")
# sns.kdeplot(df_lse["log10_Ks"], label="LSE", color="seagreen", fill=False, linewidth=1, alpha=0.6)
# sns.kdeplot(df_whole["log10_Ks"], label="Whole Genome", color="gray", fill=False, linewidth=1, alpha=0.6)
plt.xlabel("log₁₀(Kₛ)", fontsize=14)
plt.ylabel("Gene-pair Frequency", fontsize=14)
plt.ylim(0,2)
# plt.title("Synonymous Substitution Rate (log₁₀Kₛ)", fontsize=15)
plt.legend(title="Category",loc="upper right",frameon=False)
plt.tight_layout()
plt.savefig("overlay_log10Ks_lse_vs_genome.png", dpi=300)




# --- Plotting: log10(omega) overlay ---
plt.figure(figsize=(6, 5))
bins = np.linspace(-2.99, 1, 30)
plt.hist(df_lse["log10_omega"], bins=bins, alpha=0.6, label="LSE", density=True, color="seagreen")
plt.hist(df_whole["log10_omega"], bins=bins, alpha=0.6, label="Whole Genome", density=True,color="gray")
# sns.kdeplot(df_lse["log10_omega"], label="LSE", color="seagreen", fill=False, linewidth=1, alpha=0.6)
# sns.kdeplot(df_whole["log10_omega"], label="Whole Genome", color="gray", fill=False, linewidth=1, alpha=0.6)
plt.xlabel("log₁₀(ω = Kₐ/Kₛ)", fontsize=14)
plt.ylabel("Gene-pair Frequency", fontsize=14)
plt.ylim(0,2)
# plt.title("Selection Pressure (log₁₀ω)", fontsize=15)
plt.legend(title="Category",loc="upper right",frameon=False)
plt.tight_layout()
plt.savefig("overlay_log10omega_lse_vs_genome.png", dpi=300)



# --- Plotting: log10(Ka) overlay ---
plt.figure(figsize=(6, 5))
bins = np.linspace(-2.99,1, 30)
plt.hist(df_lse["log10_Ka"], bins=bins, alpha=0.6, label="LSE", density=True, color="seagreen")
plt.hist(df_whole["log10_Ka"], bins=bins, alpha=0.6, label="Whole Genome", density=True, color="gray")
# sns.kdeplot(df_lse["log10_Ka"], label="LSE", color="seagreen", fill=False, linewidth=1, alpha=0.6)
# sns.kdeplot(df_whole["log10_Ka"], label="Whole Genome", color="gray", fill=False, linewidth=1, alpha=0.6)
plt.xlabel("log₁₀(Kₐ)", fontsize=14)
plt.ylabel("Gene-pair Frequency", fontsize=14)
plt.ylim(0,2)
# plt.title("Nonsynonymous Substitution Rate (log₁₀Kₐ)", fontsize=15)
plt.legend(title="Category",loc="upper right",frameon=False)
plt.tight_layout()
plt.savefig("overlay_log10Ka_lse_vs_genome.png", dpi=300)



from scipy.stats import mannwhitneyu, ks_2samp

def report_tests(feature, label=""):
    print(f"\n=== {label} ({feature}) ===")

    # Mann–Whitney U test (tests for shift in central tendency)
    u_stat, u_p = mannwhitneyu(df_lse[feature], df_whole[feature], alternative="two-sided")
    print(f"Mann–Whitney U test: U = {u_stat:.2e}, p = {u_p:.3e}")

    # Kolmogorov–Smirnov test (tests for difference in distributions)
    ks_stat, ks_p = ks_2samp(df_lse[feature], df_whole[feature])
    print(f"Kolmogorov–Smirnov test: D = {ks_stat:.3f}, p = {ks_p:.3e}")

# Run tests for each indicator
report_tests("omega", label="Selection Pressure")
report_tests("Ks", label="Synonymous Substitution Rate")
report_tests("Ka", label="Nonsynonymous Substitution Rate")

def cliffs_delta(x, y):
    """
    Compute Cliff's Delta (δ) between two samples.
    Returns delta and the interpreted magnitude.
    """
    from itertools import product

    n_x, n_y = len(x), len(y)
    more = sum(1 for xi, yi in product(x, y) if xi > yi)
    less = sum(1 for xi, yi in product(x, y) if xi < yi)
    delta = (more - less) / (n_x * n_y)

    # Interpret the magnitude
    abs_delta = abs(delta)
    if abs_delta < 0.147:
        magnitude = "negligible"
    elif abs_delta < 0.33:
        magnitude = "small"
    elif abs_delta < 0.474:
        magnitude = "medium"
    else:
        magnitude = "large"

    return delta, magnitude

delta_ks, mag_ks = cliffs_delta(df_lse["Ks"], df_whole["Ks"])
delta_ka, mag_ka = cliffs_delta(df_lse["Ka"], df_whole["Ka"])
delta_omega, mag_omega = cliffs_delta(df_lse["omega"], df_whole["omega"])

print(f"Ks: Cliff’s delta = {delta_ks:.3f} ({mag_ks})")
print(f"Ka: Cliff’s delta = {delta_ka:.3f} ({mag_ka})")
print(f"ω:  Cliff’s delta = {delta_omega:.3f} ({mag_omega})")


# --- Load WGD Ks result ---
df = pd.read_csv("lse_900genes.fa.tsv.ks.tsv", sep="\t")
df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=["dN", "dS", "dN/dS", "family", "g1", "g2"])
df = df[(df["dN"] > 0) & (df["dS"] > 0) & (df["dN/dS"] > 0) & (df["dN/dS"] < 10)]

# --- Load BED file with gene names ---
bed_cols = ["chrom", "start", "end", "gene_id", "gene_name","strand"]
gene_info = pd.read_csv("genes_function.bed", sep="\t", names=bed_cols)

# Strip ".t1" suffix if present
gene_info["gene_id"] = gene_info["gene_id"].str.replace(r"\.t\d+$", "", regex=True)

# Build a mapping from gene_id to gene_name
gene_id_to_name = dict(zip(gene_info["gene_id"], gene_info["gene_name"]))

# Re-map gene names (after fixing g1/g2 formatting as shown earlier)
# Strip off .t1 and split the pair into two gene IDs
df[["g1", "g2"]] = df["pair"].str.replace(r"\.t\d+$", "", regex=True).str.split("__", expand=True)
# If your BED gene_info DataFrame has these IDs in "gene_id"
gene_id_to_name = dict(zip(gene_info["gene_id"], gene_info["gene_name"]))

# Map names to each gene
df["gene1_name"] = df["g1"].map(gene_id_to_name)
df["gene2_name"] = df["g2"].map(gene_id_to_name)


# Create list of unique gene names per family
def collect_gene_names(subdf):
    names = pd.concat([subdf["gene1_name"], subdf["gene2_name"]]).dropna().unique()
    return ", ".join(sorted(names))

# --- Summary: omega stats ---
summary_stats = (
    df.groupby("family")
    .agg(
        mean_omega=("dN/dS", "mean"),
        median_omega=("dN/dS", "median"),
        max_omega=("dN/dS", "max"),
        n_pairs=("dN/dS", "count")
    )
    .reset_index()
)

# --- Compute representative gene names separately ---
def collect_gene_names(subdf):
    names = pd.concat([subdf["gene1_name"], subdf["gene2_name"]]).dropna().unique()
    return ", ".join(sorted(names))

gene_name_map = df.groupby("family").apply(collect_gene_names).reset_index()
print(gene_name_map)
gene_name_map.columns = ["family", "representative_genes"]

# --- Merge and sort ---
summary = pd.merge(summary_stats, gene_name_map, on="family")
summary = summary.sort_values(by="mean_omega", ascending=False)

# --- Save ---
summary.to_csv("LSE_gene_family_omega_summary_with_names.tsv", sep="\t", index=False)
# print(summary.head(10))


















