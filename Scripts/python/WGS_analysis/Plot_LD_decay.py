import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/snp_analysis")

# === Load LD summary ===
df = pd.read_csv("ld_decay_summary.tsv", sep="\t")
# === Filter low-confidence bins ===
min_pair_count = 100
df = df[df["Count"] >= min_pair_count]

# === Set up plot ===
fig, ax1 = plt.subplots(figsize=(6,5))

ax1.plot(df["Distance"], df["Median_r2"],
         color="darkblue", lw=2, label="Median $r^2$")

# -------- set limits here ----------
ax1.set_xlim(0, 250_000)      # show only 0-250 kb on the x-axis
ax1.set_ylim(0, 0.6)        # y-axis from 0.20 to 0.60
# -----------------------------------

ax1.set_xlabel("Distance (bp)")
ax1.set_ylabel("Median $r^2$", color="darkblue")
ax1.grid(True, linestyle="--", alpha=0.5)

ax2 = ax1.twinx()
ax2.fill_between(df["Distance"], df["Count"],
                 color="lightgray", alpha=0.5, step="mid")

ax2.set_ylabel("Pair count", color="gray", fontsize=12)
ax2.tick_params(axis="y", labelcolor="gray")

# Titles and layout
# plt.title("LD decay in C.virginica", fontsize=14)
fig.tight_layout()

# Save and show
plt.savefig("ld_decay_plot.png", dpi=300)
# plt.show()



recomb_rate_per_bp = 3e-8         # 1 cM/Mb (adjust if known)
# === Estimate effective population size ===
df["c"] = df["Distance"] * recomb_rate_per_bp

# Avoid divide-by-zero
df = df[df["c"] > 0]

# Estimate Ne from r²
df["Ne_median"] = (1 / (4 * df["c"])) * (1 / df["Median_r2"] - 1)
df["Ne_mean"] = (1 / (4 * df["c"])) * (1 / df["Mean_r2"] - 1)
print(df["Ne_median"])
# === Plot ===
plt.figure(figsize=(6, 5))
plt.plot(df["Distance"], df["Ne_median"], label="Ne (median r²)", color="blue")
# plt.plot(df["Distance"], df["Ne_mean"], label="Ne (mean r²)", color="green", alpha=0.6)
plt.xlabel("Distance between SNPs (bp)")
plt.ylabel("Estimated effective population size (Ne)")
# plt.title("Estimated Ne from LD decay")
# plt.legend()
plt.ylim(0, df[["Ne_median", "Ne_mean"]].quantile(0.95).max())  # Optional to reduce y-axis outliers
plt.grid(True)
plt.tight_layout()
plt.savefig("Ne_from_LD_decay.png", dpi=300)
# plt.show()

# === Save results ===
df.to_csv("ld_decay_with_Ne.tsv", sep="\t", index=False)


# Filter for analysis window (e.g., first 250kb)
ld_summary = df[df["Distance"] <= 250_000]

# Calculate basic statistics
summary_stats = {
    "Mean rsquare": ld_summary["Median_r2"].mean(),
    "Median rsquare": ld_summary["Median_r2"].median(),
    "Min rsquare": ld_summary["Median_r2"].min(),
    "Max rsquare": ld_summary["Median_r2"].max(),
    "Distance at which median rsquare < 0.2": ld_summary[ld_summary["Median_r2"] < 0.2]["Distance"].min(),
    "Mean pair count per bin": ld_summary["Count"].mean(),
    "Total SNP pairs analyzed": ld_summary["Count"].sum(),
    "Mean estimated Ne (median rsquare)": ld_summary["Ne_median"].mean(),
    "Median estimated Ne (median rsquare)": ld_summary["Ne_median"].median()
}

# Format as DataFrame
import pandas as pd
panel3_summary = pd.DataFrame.from_dict(summary_stats, orient="index", columns=["Value"])
panel3_summary.index.name = "Statistic"

# Save to file
panel3_summary.to_csv("panel3_LD_summary.tsv", sep="\t")

# Print for review
print(panel3_summary)