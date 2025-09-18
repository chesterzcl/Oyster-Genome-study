#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import os
import numpy as np

# ===============================
# Setup & I/O
# ===============================
working_dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2"
os.chdir(working_dir)

# ===============================
# 1. Load TF class/family mapping
# ===============================
parsed_entries = []
with open("TFclass.txt") as f:
    entry = {}
    for line in f:
        line = line.strip()
        if line.startswith("AC "):
            entry['matrix_id'] = line.split()[1]
        elif line.startswith("ID "):
            entry['TF_name'] = line.split()[1]
        elif line.startswith("CC tf_family:"):
            entry['TF_family'] = line.split(":", 1)[1].strip() or "Unknown"
        elif line.startswith("CC tf_class:"):
            entry['TF_class'] = line.split(":", 1)[1].strip() or "Unknown"
        elif line == "//":
            if 'matrix_id' in entry:
                parsed_entries.append(entry)
            entry = {}

mapping_df = pd.DataFrame(parsed_entries).fillna("Unknown")
print("[INFO] Loaded mapping:", mapping_df.shape)
mapping_df.to_csv("motif_family_mapping.csv", index=False)


# ===============================
# 2. Load scored motif calls
# ===============================
cols = ["chr", "start", "end", "motif_name", "pwm_score", "strand", "footprint_score"]
df = pd.read_csv("ScoredMotifs.bed", sep="\t", header=None, names=cols)
print("[INFO] Loaded motif data:", df.shape)


# ===============================
# 3. QC and filtering
# ===============================
nonzero = df[df["footprint_score"] > 0]
threshold = 1.5
high_conf = nonzero[nonzero["footprint_score"] >= threshold].copy()
print("[INFO] High-confidence calls:", high_conf.shape)
print(nonzero["footprint_score"].describe())
nonzero["log_footprint_score"] = np.log1p(nonzero["footprint_score"])
plt.figure(figsize=(8,5))
plt.hist(nonzero["log_footprint_score"], bins=100, color='gray', edgecolor='black', alpha=0.7)
plt.xlabel("log(1 + Footprint Score)")
plt.ylabel("Count")
plt.title("Log-transformed Footprint Score Distribution")
plt.tight_layout()
plt.show()


# Add matrix_id for merge
high_conf['matrix_id']=high_conf['motif_name'].apply(lambda x: x.split('_')[-1])
high_conf_annot=high_conf.merge(mapping_df,on='matrix_id',how='left')
high_conf_annot['TF_family']=high_conf_annot['TF_family'].fillna("Unknown")
high_conf_annot['TF_class']=high_conf_annot['TF_class'].fillna("Unknown")
high_conf_annot.to_csv("SE_footprint_raw.bed", sep="\t", header=True,index=False)
print("[INFO] Annotated calls:",high_conf_annot.shape)
print(list(high_conf_annot))
print(high_conf_annot)
print(set(high_conf_annot['TF_class']))
print(len(set(high_conf_annot['TF_class'])))


def plot_tf_class_signal_from_unmerged(
    high_conf_annot,
    region_chr,
    region_start,
    region_end,
    bin_size=50,
    aggfunc='sum',
    min_total_score=0,
    save_path=None,
    figsize=(12, 5),
    class_palette='tab20'
):
    # Filter region
    region_data = high_conf_annot[
        (high_conf_annot['chr'] == region_chr) &
        (high_conf_annot['start'] <= region_end) &
        (high_conf_annot['end'] >= region_start)
    ].copy()

    if region_data.empty:
        print(f"[WARN] No TF footprints in region {region_chr}:{region_start}-{region_end}")
        return None

    print(f"[INFO] Found {len(region_data)} unmerged TF hits in region.")

    # Define bins
    bins = np.arange(region_start, region_end + bin_size, bin_size)
    bin_centers = bins[:-1] + bin_size // 2

    classes = sorted(region_data['TF_class'].dropna().unique())
    all_records = []

    # Aggregate signal per bin per class
    for tf_class in classes:
        class_hits = region_data[region_data['TF_class'] == tf_class]
        scores_per_bin = []
        for b_start, b_end in zip(bins[:-1], bins[1:]):
            overlaps = class_hits[
                (class_hits['start'] < b_end) & (class_hits['end'] > b_start)
            ]
            if aggfunc == 'sum':
                score = overlaps['footprint_score'].sum()
            else:
                score = overlaps['footprint_score'].mean() if not overlaps.empty else 0
            scores_per_bin.append(score)

        total_score = sum(scores_per_bin)
        if total_score >= min_total_score:
            for x, y in zip(bin_centers, scores_per_bin):
                all_records.append({
                    'Bin': x,
                    'TF_class': tf_class,
                    'Footprint_Signal': y
                })

    if not all_records:
        print("[WARN] No signal to plot after aggregation/filtering.")
        return None

    # Build DataFrame
    plot_df = pd.DataFrame(all_records)

    # Plotting
    plt.figure(figsize=figsize)
    ax = sns.lineplot(
        data=plot_df,
        x='Bin',
        y='Footprint_Signal',
        hue='TF_class',
        palette=class_palette
    )

    ax.set_xlabel(f"Genomic coordinate ({region_chr})", fontsize=12)
    ax.set_ylabel("Footprint signal", fontsize=12)
    # ax.set_title(f"TF Class Footprint Signal in {region_chr}:{region_start}-{region_end}",fontsize=14)
    # ax.legend(title='TF Class', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8)
    sns.despine()
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"[INFO] Saved track plot to {save_path}")
    else:
        plt.show()

    return plot_df


plot_tf_class_signal_from_unmerged(
    high_conf_annot,
    region_chr='HiC_scaffold_3',
    region_start=37477000,
    region_end=37477400,
    bin_size=50,
    aggfunc='mean',
    min_total_score=0,
    save_path='TF_Footprint_Class_Signal_Track_Unmerged.png'
)


# ===============================
# 4. Define merging functions
# ===============================
def merge_family_hits_with_class_and_scores(df, family_class_df, merge_distance=10):
    merged_all = []
    df['TF_family'] = df['TF_family'].fillna("Unknown")

    for family, group in df.groupby('TF_family'):
        if (family_class_df['TF_family'] == family).any():
            tf_class = family_class_df.query('TF_family == @family')['TF_class'].iloc[0]
        else:
            tf_class = 'Unknown'
        group = group.sort_values(['chr', 'start'])

        current_chr = None
        current_start = None
        current_end = None
        pwm_scores = []
        footprint_scores = []

        for _, row in group.iterrows():
            chrom, start, end = row['chr'], row['start'], row['end']
            pwm, footprint = row['pwm_score'], row['footprint_score']

            if current_chr != chrom:
                if current_chr is not None:
                    merged_all.append([current_chr, current_start, current_end, family, tf_class,
                                       np.max(pwm_scores), np.mean(footprint_scores)])
                current_chr, current_start, current_end = chrom, start, end
                pwm_scores, footprint_scores = [pwm], [footprint]
            elif start <= current_end + merge_distance:
                current_end = max(current_end, end)
                pwm_scores.append(pwm)
                footprint_scores.append(footprint)
            else:
                merged_all.append([current_chr, current_start, current_end, family, tf_class,
                                   np.max(pwm_scores), np.mean(footprint_scores)])
                current_chr, current_start, current_end = chrom, start, end
                pwm_scores, footprint_scores = [pwm], [footprint]

        if current_chr is not None:
            merged_all.append([current_chr, current_start, current_end, family, tf_class,
                               np.max(pwm_scores), np.mean(footprint_scores)])

    merged_df = pd.DataFrame(merged_all, columns=[
        'chr', 'start', 'end', 'TF_family', 'TF_class', 'max_pwm_score', 'mean_footprint_score'
    ])
    return merged_df.sort_values(['chr', 'start']).reset_index(drop=True)


def pick_best_interval(group):
    top_footprint = group['mean_footprint_score'].max()
    best = group[group['mean_footprint_score'] == top_footprint]
    if len(best) > 1:
        top_pwm = best['max_pwm_score'].max()
        best = best[best['max_pwm_score'] == top_pwm]
    return best.iloc[0]


def resolve_cross_family_overlaps(merged_df):
    final_intervals = []
    for chrom in merged_df['chr'].unique():
        chrom_df = merged_df[merged_df['chr'] == chrom].sort_values('start').reset_index(drop=True)
        active_group = []
        last_end = None
        for _, row in chrom_df.iterrows():
            if not active_group:
                active_group = [row]
                last_end = row['end']
            elif row['start'] <= last_end:
                active_group.append(row)
                last_end = max(last_end, row['end'])
            else:
                best = pick_best_interval(pd.DataFrame(active_group))
                final_intervals.append(best)
                active_group = [row]
                last_end = row['end']
        if active_group:
            best = pick_best_interval(pd.DataFrame(active_group))
            final_intervals.append(best)

    return pd.DataFrame(final_intervals).sort_values(['chr', 'start']).reset_index(drop=True)


# ===============================
# 5. Merging within families
# ===============================
print("\n[Merging] Step 1: Collapse hits within families ...")
merged_family_sites = merge_family_hits_with_class_and_scores(high_conf_annot, mapping_df, merge_distance=10)
merged_family_sites.to_csv("SE_Footprints_MergedFamilySites_WithScoresAndClass.bed", sep="\t", index=False)
print("[INFO] Merged family-level intervals:", merged_family_sites.shape)


# ===============================
# 6. Removing residual overlaps
# ===============================
print("\n[De-overlapping] Step 2: Removing residual overlaps across families ...")
final_df = resolve_cross_family_overlaps(merged_family_sites)
final_df.to_csv("SE_Footprints_Final_NonOverlapping_Sites.bed", sep="\t", index=False)
print("[INFO] Final non-overlapping intervals:", final_df.shape)


# ===============================
# 7. Summarize and plot results
# ===============================
family_summary_merged = final_df.groupby('TF_family').agg(
    n_sites=('TF_family','count'),
    mean_pwm_score=('max_pwm_score','mean'),
    mean_footprint_score=('mean_footprint_score','mean')
).sort_values('n_sites',ascending=False)
family_summary_merged.to_csv("SE_Footprint_Summary_by_TF_Family_Merged.csv")


class_summary_merged = final_df.groupby('TF_class').agg(
    n_sites=('TF_class', 'count'),
    mean_pwm_score=('max_pwm_score', 'mean'),
    mean_footprint_score=('mean_footprint_score', 'mean')
).sort_values('n_sites', ascending=False)
class_summary_merged.to_csv("SE_Footprint_Summary_by_TF_Class_Merged.csv")


top_n=20
plt.figure(figsize=(10,6))
sns.barplot(data=family_summary_merged.head(top_n).reset_index(), y="TF_family", x="n_sites", color='forestgreen')
plt.title(f"Top {top_n} TF families by *unique* merged binding sites (score ≥ {threshold})")
plt.xlabel("Number of unique binding sites")
plt.ylabel("TF Family")
plt.tight_layout()
plt.savefig("Step9_TopTF_Families_ByCount_Merged.png")
plt.close()


plt.figure(figsize=(10,6))
sns.barplot(data=class_summary_merged.head(top_n).reset_index(), y="TF_class", x="n_sites", color='slateblue')
plt.title(f"Top {top_n} TF classes by *unique* merged binding sites (score ≥ {threshold})")
plt.xlabel("Number of unique binding sites")
plt.ylabel("TF Class")
plt.tight_layout()
plt.savefig("Step10_TopTF_Classes_ByCount_Merged.png")
plt.close()


# Output unique families
unique_families = final_df['TF_family'].dropna().unique()
unique_families.sort()
pd.DataFrame({'TF_family': unique_families}).to_csv("Unique_TF_Families_Sorted.csv", index=False)


print("\n[Done!] All outputs generated successfully.")