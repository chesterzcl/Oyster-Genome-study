import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster")
# Load eigenvectors (output from PLINK)
evec = pd.read_csv("cv20_maf05_pca.eigenvec", delim_whitespace=True, header=None)

# Assign column names (first two columns are FID and IID, rest are PC1-PCn)
evec.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, evec.shape[1] - 1)]

# Extract short ID (e.g., "10L") and group ("Large" or "Small")
def extract_short_id(iid):
    match = re.match(r"(\d+[LS])_", iid)
    return match.group(1) if match else iid

evec["ShortID"] = evec["IID"].apply(extract_short_id)
evec["Group"] = evec["ShortID"].apply(lambda x: "Large" if x.endswith("L") else "Small")

# Plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=evec, x="PC1", y="PC2", hue="Group", s=80, edgecolor="black")

# Annotate each point with ShortID
for _, row in evec.iterrows():
    plt.text(row["PC1"], row["PC2"] + 0.012, row["ShortID"], fontsize=6, ha="center")

plt.title("PCA for first 20 samples")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(title="Group", loc="upper right")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("oyster_pca_plot_grouped.png", dpi=300)
plt.show()