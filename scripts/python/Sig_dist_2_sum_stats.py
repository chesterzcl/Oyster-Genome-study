import pandas as pd
import re
from scipy.stats import fisher_exact
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster")

ip_file = "20cv_df4_hoh_biSNP_filtered_ann_sig_dist_gene.txt"

df = pd.read_csv(ip_file, sep='\t', header=0)

# === Extract sample columns ===
sample_cols = df.columns[5:]  # First 5 are assumed: Chromosome, position, ref, alt, annotation
l_columns = [col for col in sample_cols if 'L' in col]
s_columns = [col for col in sample_cols if 'S' in col]

# === Genotype parsing function ===
def parse_genotype(entry):
    if pd.isna(entry) or entry.startswith('./.') or entry == '.' or entry == './.':
        return (0, 0)
    try:
        gt = entry.split(";")[0]
        alleles = re.split(r'[\/|]', gt)
        total = len(alleles)
        alt = alleles.count('1')
        return (alt, total)
    except:
        return (0, 0)

# === Count function per group ===
def count_group_alleles(row, group_cols):
    alt_total = 0
    allele_total = 0
    for col in group_cols:
        alt, total = parse_genotype(row[col])
        alt_total += alt
        allele_total += total
    return alt_total, allele_total

# === Apply counting and Fisher's test ===
l_freq_list, s_freq_list, fisher_pvalues = [], [], []

for idx, row in df.iterrows():
    l_alt, l_total = count_group_alleles(row, l_columns)
    s_alt, s_total = count_group_alleles(row, s_columns)
    l_ref = l_total - l_alt
    s_ref = s_total - s_alt

    l_freq_list.append(f"{l_alt}/{l_total}")
    s_freq_list.append(f"{s_alt}/{s_total}")

    if (l_total > 0 and s_total > 0):
        table = [[l_alt, l_ref], [s_alt, s_ref]]
        _, p_value = fisher_exact(table)
    else:
        p_value = 1.0
    fisher_pvalues.append(p_value)

# === Add results to dataframe ===
df['Large_group_alt_allele_freq'] = l_freq_list
df['Small_group_alt_allele_freq'] = s_freq_list
df['fisher_p_val'] = fisher_pvalues

# === Save top 20 most significant ===
top_sig_df = df.sort_values(by='fisher_p_val').head(20)

print(top_sig_df)

top_sig_df.to_csv("top20_significant_snps_gene.txt", sep='\t', index=False, columns=[
    '#Chromosome', 'position', 'ref_allele', 'alt_allele', 'annotation',
    'Large_group_alt_allele_freq', 'Small_group_alt_allele_freq', 'fisher_p_val'
])

# === Save full annotated table ===
df.to_csv("allele_freq_fisher_results_gene.txt", sep='\t', index=False)