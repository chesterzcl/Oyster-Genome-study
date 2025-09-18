import pandas as pd
from scipy.stats import fisher_exact

# === Load CNV matrix ===
matrix = pd.read_csv("cnv_cluster_presence_matrix.tsv", sep="\t", index_col=0)

# === Parse group assignments from column names ===
sample_names = matrix.columns

# Determine groups based on pattern (L or S before ".SM")
group_labels = {}
for sample in sample_names:
    if "-L.SM" in sample or "-1L.SM" in sample or "-2L.SM" in sample or "-3L.SM" in sample:
        group_labels[sample] = "L"
    elif "-S.SM" in sample:
        group_labels[sample] = "S"
    else:
        # fallback: parse between '-' and '.SM'
        token = sample.split('-')[-1].split('.')[0]
        group_labels[sample] = "L" if "L" in token else "S"

# Convert to pandas Series
groups = pd.Series(group_labels)

# Filter matrix to only include samples with group info
matrix = matrix.loc[:, matrix.columns.isin(groups.index)]
groups = groups.loc[matrix.columns]

# Define sample lists
group_S = groups[groups == "S"].index
group_L = groups[groups == "L"].index

# === Fisher's exact test per CNV cluster ===
results = []

for cluster_id, row in matrix.iterrows():
    s_present = row[group_S].sum()
    s_absent = len(group_S) - s_present
    l_present = row[group_L].sum()
    l_absent = len(group_L) - l_present

    table = [[s_present, s_absent], [l_present, l_absent]]
    
    try:
        oddsratio, pvalue = fisher_exact(table, alternative='two-sided')
    except:
        oddsratio, pvalue = float('nan'), 1.0

    results.append({
        "cluster_id": cluster_id,
        "S_present": s_present,
        "S_absent": s_absent,
        "L_present": l_present,
        "L_absent": l_absent,
        "odds_ratio": oddsratio,
        "p_value": pvalue
    })

# === Output results ===
# === Prepare BED-formatted output ===
import re

bed_records = []

for row in results:
    cluster_id = row["cluster_id"]
    match = re.match(r"(HiC_scaffold_\d+):(\d+)-(\d+)", cluster_id)
    if not match:
        continue
    chrom, start, end = match.groups()
    start, end = int(start), int(end)
    name = cluster_id
    score = f'{row["p_value"]:.4g}'  # use p-value as BED "score" column

    bed_records.append([
        chrom, start, end, name,
        row["S_present"], row["S_absent"],
        row["L_present"], row["L_absent"],
        f'{row["odds_ratio"]:.4g}', f'{row["p_value"]:.4g}'
    ])

# Convert to DataFrame
bed_df = pd.DataFrame(bed_records, columns=[
    "chrom", "start", "end", "cluster_id",
    "S_present", "S_absent", "L_present", "L_absent",
    "odds_ratio", "p_value"
])

# Sort by chromosome number and start
bed_df["chrom_num"] = bed_df["chrom"].str.extract(r'_(\d+)').astype(int)
bed_df.sort_values(["chrom_num", "start"], inplace=True)
bed_df.drop(columns="chrom_num", inplace=True)

# Output as BED-style file
bed_df.to_csv("cnv_fisher_test_S_vs_L.bed", sep="\t", index=False, header=False)
