import os
import re
from collections import Counter
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Set working directory
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/Orthogroups/")

# Load the data (skip comment lines)
df = pd.read_csv("Base_clade_results.txt", sep="\t")

# Extract species/clade names from Taxon_ID
df["Label"] = df["Taxon_ID"].str.extract(r'([A-Za-z\.]+)')[0]

# Drop rows where label could not be extracted (e.g., <6>, <7>)
df = df.dropna(subset=["Label"])

# Make sure labels are strings
df["Label"] = df["Label"].astype(str)

# Optional: define label order for sorting
label_order = ["M.angulata", "M.gigas",  "C.virginica","O.edulis"]
df["Label"] = pd.Categorical(df["Label"], categories=label_order, ordered=True)
df = df.sort_values("Label")


plt.figure(figsize=(8, 5))
plt.barh(df["Label"], df["Increase"], label="Expansion", color="steelblue")
plt.barh(df["Label"], df["Decrease"], left=df["Increase"], label="Contraction", color="indianred")

# Labeling
plt.xlabel("Number of Gene Families")
# plt.title("Gene Family Expansions and Contractions per Species (CAFE5)")
plt.legend(title="", loc="upper right",frameon=False)
sns.despine(top=True,right=True,left=True)
plt.tight_layout()

# Save and show
plt.savefig("cafe5_expansion_contraction_horizontal.png", dpi=300)
plt.show()


