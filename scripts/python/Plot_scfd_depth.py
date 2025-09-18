import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.signal import find_peaks


INPUT_FILE="/Users/lizicheng/Desktop/Data/dog_genome/oyster/depth_distribution_by_scaffold_wgs.tsv" 
OUTPUT_DIR="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/coverage_plot/"
OUTPUT_FILE =OUTPUT_DIR+"/coverage_summary_wgs.tsv"
DEPTH_THRESHOLD =500    # Show coverage distribution for x < 500 only
# -------------------------------


os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load file
df=pd.read_csv(INPUT_FILE, sep="\t")
df["Depth"]=pd.to_numeric(df["Depth"], errors="coerce")
df["Count"]=pd.to_numeric(df["Count"], errors="coerce")
df.dropna(inplace=True)


# Weighted percentile function
def weighted_percentile(data, percentiles):
    data_sorted = data.sort_values("Depth")
    cumsum = data_sorted["Count"].cumsum()
    total = cumsum.iloc[-1]
    results = []
    for p in percentiles:
        cutoff = total * p / 100.0
        depth = data_sorted.loc[cumsum >= cutoff, "Depth"].iloc[0]
        results.append(depth)
    return results

# Weighted standard deviation
def weighted_std(depths, weights, mean):
    return np.sqrt(np.average((depths - mean) ** 2, weights=weights))

# Per-scaffold summary
def summarize(group):
    total = group["Count"].sum()
    mean = (group["Depth"] * group["Count"]).sum() / total
    stddev = weighted_std(group["Depth"], group["Count"],mean)
    cv = stddev / mean if mean > 0 else 0
    median, p10, p25, p75, p90 = weighted_percentile(group, [50, 10, 25, 75, 90])
    return pd.Series({
        "Mean_depth": round(mean, 1),
        "Std_depth":round(stddev,1),
        "CV": round(cv, 3),
        "Median_Depth": round(median, 1),
        "10th_Percentile": round(p10, 1),
        "25th_Percentile": round(p25, 1),
        "75th_Percentile": round(p75, 1),
        "90th_Percentile": round(p90, 1)
    })

summary = df.groupby("Scaffold").apply(summarize).reset_index()

# Genome-wide summary
genome_row = summarize(df)
genome_row["Scaffold"] = "Genome-wide"
summary = pd.concat([summary, genome_row.to_frame().T], ignore_index=True)

# Save output
summary.to_csv(OUTPUT_FILE, sep="\t", index=False)


# Filter by depth
df=df[df["Depth"] < DEPTH_THRESHOLD]

for scaffold in df["Scaffold"].unique():
    sub = df[df["Scaffold"] == scaffold].sort_values("Depth")
    total = sub["Count"].sum()
    sub["Frequency"] = sub["Count"] / total

    scaffold_id=scaffold.split('_')[2]

    plt.figure(figsize=(10, 4))
    plt.plot(sub["Depth"], sub["Frequency"],color="darkgreen",linewidth=1.5)
    plt.title(f"Depth Frequency Distribution:Chr{scaffold_id}")
    plt.xlabel("Depth")
    plt.ylabel("Frequency (Proportion of Bases)")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/{scaffold}_depth_freq_wgs.png")
    plt.close()


