import pandas as pd
import numpy as np
from collections import Counter
from scipy.stats import chi2_contingency, fisher_exact

dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/"
vcf_path = "merged_sv.sorted.vcf"

group_L = ["10L.SM", "11L.SM", "12L.SM", "3L.SM", "4L.SM", "5L.SM", "6L.SM", "7L.SM", "8L.SM", "9L.SM"]
group_S = ["11S.SM", "12S.SM", "14S.SM", "1S.SM", "2S.SM", "3S.SM", "4S.SM", "6S.SM", "8S.SM", "9S.SM"]

min_len, max_len = 50, 1_000_000
min_supp = 16
valid_svtypes = {"DEL", "DUP", "INS", "INV"}

sample_ids = []
results = []

with open(dir+vcf_path) as f:
    for line in f:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            header_fields = line.strip().split("\t")
            sample_ids = header_fields[9:]
            sample_index = {s: i for i, s in enumerate(sample_ids)}
            continue

        fields = line.strip().split("\t")
        chrom, pos, _, ref, alt, _, _, info_str, format_str = fields[:9]
        samples = fields[9:]
        pos = int(pos)

        # Parse INFO fields
        info = dict(item.split("=", 1) if "=" in item else (item, True) for item in info_str.split(";"))
        svtype = info.get("SVTYPE", "")
        if svtype not in valid_svtypes:
            continue
        try:
            end = int(info.get("END", pos))
            svlen = abs(int(info.get("SVLEN", end - pos)))
            supp = int(info.get("SUPP", 0))
        except ValueError:
            continue
        if not (min_len < svlen < max_len) or supp < min_supp:
            continue

        # Parse FORMAT
        format_keys = format_str.split(":")
        if "GT" not in format_keys:
            continue
        gt_index = format_keys.index("GT")

        def extract_gt(sample_str):
            try:
                return sample_str.split(":")[gt_index]
            except:
                return "./."

        L_gts = [extract_gt(samples[sample_index[s]]) for s in group_L if s in sample_index]
        S_gts = [extract_gt(samples[sample_index[s]]) for s in group_S if s in sample_index]

        def count_genotypes(gt_list):
            return Counter(gt for gt in gt_list if gt != "./.")

        L_ct = count_genotypes(L_gts)
        S_ct = count_genotypes(S_gts)

        all_gts = sorted(set(L_ct.keys()) | set(S_ct.keys()))
        if len(all_gts) < 2:
            continue

        table = [
            [L_ct.get(gt, 0) for gt in all_gts],
            [S_ct.get(gt, 0) for gt in all_gts]
        ]

        # Perform test
        try:
            if all(v < 5 for row in table for v in row) and len(all_gts) == 2:
                _, pval = fisher_exact(table)
            else:
                _, pval, _, _ = chi2_contingency(table)
        except:
            pval = np.nan

        if pval < 0.05:
            results.append({
                "CHROM": chrom,
                "POS": pos,
                "END": end,
                "SVTYPE": svtype,
                "SVLEN": svlen,
                "SUPP": supp,
                "L_counts": dict(L_ct),
                "S_counts": dict(S_ct),
                "PVAL": pval,
                "-log10(PVAL)": -np.log10(pval)
            })

# Save only significant ones
sig_df = pd.DataFrame(results)
sig_df.to_csv(dir+"sv_genotype_comparison.txt", sep="\t", index=False)

print(f"{len(sig_df)} SVs tested. {sum(sig_df['PVAL'] < 0.05)} had p < 0.05.")