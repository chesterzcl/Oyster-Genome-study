import pybedtools
from pybedtools import BedTool
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
from matplotlib.lines import Line2D

# === Input files ===
dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/bed/"
retrogene_bed = "retrogenes.bed"
nonretrogene_bed = "non_retrogenes.bed"
retrotransposon_bed = "retroelements.bed"
atac_high_file = "ATAC_high_variance_peaks.bed"
atac_low_file = "ATAC_low_variance_peaks.bed"
proximity=100

prox_list=list(range(1,1251,50))


rg_te_proximity_vec=[]
nrg_te_proximity_vec=[]
P_val_list=[]

for proximity in prox_list:

	retro_df = pd.read_csv(dir+retrotransposon_bed, sep="\t", header=None,
	                       names=["chrom", "start", "end", "name", "dot", "class_subclass"])

	# Extract class only (e.g., LINE from LINE/RTE-X)
	retro_df["class"] = retro_df["class_subclass"].str.split("/").str[0]
	# print(retro_df)
	# Export BEDs per class
	for cls in ["LINE", "SINE", "LTR"]:
	    class_df = retro_df[retro_df["class"] == cls]
	    class_df[["chrom", "start", "end", "name", "dot", "class_subclass"]].to_csv(
	        dir+f"retroelements_{cls}.bed", sep="\t", index=False, header=False)

	# === Load BEDs ===
	retrogenes = BedTool(dir + retrogene_bed)
	nonretrogenes = BedTool(dir + nonretrogene_bed)
	retrotransposons = BedTool(dir + retrotransposon_bed)
	# retrotransposons = BedTool(dir + "retroelements_LTR.bed")
	atac_high = BedTool(dir + atac_high_file)
	atac_low = BedTool(dir + atac_low_file)
	atac_all = atac_high.cat(atac_low, postmerge=False)
	# print(retrogenes)
	# === De-duplication helper ===
	def bedtool_to_unique_genes(windowed_bedtool):
	    df = windowed_bedtool.to_dataframe(names=[
	        "chrom", "start", "end", "name", "score", "strand",
	        "feat_chrom", "feat_start", "feat_end", "feat_name", "feat_score", "feat_strand"
	    ])
	    df_unique = df.drop_duplicates(subset=["chrom", "start", "end", "name", "strand"])
	    bed_df = df_unique[["chrom", "start", "end", "name", "score", "strand"]]
	    return BedTool.from_dataframe(bed_df)

	# === Retrogene Proximity ===
	rgenes_with_TE_overlap = retrogenes.intersect(retrotransposons, u=True)
	rgenes_near_TE = bedtool_to_unique_genes(retrogenes.window(retrotransposons, w=proximity))
	rgenes_near_ATAC_high = bedtool_to_unique_genes(retrogenes.window(atac_high, w=proximity))
	rgenes_near_ATAC_low = bedtool_to_unique_genes(retrogenes.window(atac_low, w=proximity))
	rgenes_near_ATAC_all = bedtool_to_unique_genes(retrogenes.window(atac_all, w=proximity))

	# === Non-retrogene Proximity ===
	nongenes_with_TE_overlap = nonretrogenes.intersect(retrotransposons, u=True)
	nongenes_near_TE = bedtool_to_unique_genes(nonretrogenes.window(retrotransposons, w=proximity))
	nongenes_near_ATAC_high = bedtool_to_unique_genes(nonretrogenes.window(atac_high, w=proximity))
	nongenes_near_ATAC_low = bedtool_to_unique_genes(nonretrogenes.window(atac_low, w=proximity))
	nongenes_near_ATAC_all = bedtool_to_unique_genes(nonretrogenes.window(atac_all, w=proximity))

	# === Summary ===
	total_retrogenes = len(retrogenes)
	total_nonretrogenes = len(nonretrogenes)

	print("=== Retrogene Proximity Summary ===")
	print(f"Total retrogenes: {total_retrogenes}")
	print(f"Strict overlap with TE: {len(rgenes_with_TE_overlap)} ({len(rgenes_with_TE_overlap)/total_retrogenes:.2%})")
	print(f"Within {proximity}bp of TE: {len(rgenes_near_TE)} ({len(rgenes_near_TE)/total_retrogenes:.2%})")
	print(f"Within {proximity}bp of high-variance ATAC: {len(rgenes_near_ATAC_high)} ({len(rgenes_near_ATAC_high)/total_retrogenes:.2%})")
	print(f"Within {proximity}bp of low-variance ATAC: {len(rgenes_near_ATAC_low)} ({len(rgenes_near_ATAC_low)/total_retrogenes:.2%})")
	print(f"Within {proximity}bp of any ATAC (combined): {len(rgenes_near_ATAC_all)} ({len(rgenes_near_ATAC_all)/total_retrogenes:.2%})")

	print("\n=== Non-retrogene Proximity Summary ===")
	print(f"Total non-retrogenes: {total_nonretrogenes}")
	print(f"Strict overlap with TE: {len(nongenes_with_TE_overlap)} ({len(nongenes_with_TE_overlap)/total_nonretrogenes:.2%})")
	print(f"Within {proximity}bp of TE: {len(nongenes_near_TE)} ({len(nongenes_near_TE)/total_nonretrogenes:.2%})")
	print(f"Within {proximity}bp of high-variance ATAC: {len(nongenes_near_ATAC_high)} ({len(nongenes_near_ATAC_high)/total_nonretrogenes:.2%})")
	print(f"Within {proximity}bp of low-variance ATAC: {len(nongenes_near_ATAC_low)} ({len(nongenes_near_ATAC_low)/total_nonretrogenes:.2%})")
	print(f"Within {proximity}bp of any ATAC (combined): {len(nongenes_near_ATAC_all)} ({len(nongenes_near_ATAC_all)/total_nonretrogenes:.2%})")

	# === Optional output ===
	rgenes_with_TE_overlap.saveas(dir + "retrogenes_TE_overlap.bed")
	rgenes_near_TE.saveas(dir + "retrogenes_near_TE.bed")
	rgenes_near_ATAC_high.saveas(dir + "retrogenes_near_ATAC_high.bed")
	rgenes_near_ATAC_low.saveas(dir + "retrogenes_near_ATAC_low.bed")
	rgenes_near_ATAC_all.saveas(dir + "retrogenes_near_ATAC_combined.bed")

	nongenes_with_TE_overlap.saveas(dir + "nonretrogenes_TE_overlap.bed")
	nongenes_near_TE.saveas(dir + "nonretrogenes_near_TE.bed")
	nongenes_near_ATAC_high.saveas(dir + "nonretrogenes_near_ATAC_high.bed")
	nongenes_near_ATAC_low.saveas(dir + "nonretrogenes_near_ATAC_low.bed")
	nongenes_near_ATAC_all.saveas(dir + "nonretrogenes_near_ATAC_combined.bed")

	# === Exclude retrotransposons overlapping non-retrogenes ===
	all_genes = retrogenes.cat(nonretrogenes, postmerge=False)
	retro_free = retrotransposons.intersect(all_genes, v=True)

	# === Retrogene proximity to free TEs (within 5kb) ===
	rgenes_near_free_TE = retrogenes.window(retro_free, w=proximity)
	rgenes_near_free_TE_unique = bedtool_to_unique_genes(rgenes_near_free_TE)

	# === Non-retrogene proximity to same free TEs (within 5kb) ===
	nongenes_near_free_TE = nonretrogenes.window(retro_free, w=proximity)
	nongenes_near_free_TE_unique = bedtool_to_unique_genes(nongenes_near_free_TE)

	# === Summarize ===
	print("=== Proximity to Free (non-overlapping) TEs ===")
	print(f"Retrogenes near free TE: {len(rgenes_near_free_TE_unique)} / {len(retrogenes)} ({len(rgenes_near_free_TE_unique)/len(retrogenes):.2%})")
	print(f"Non-retrogenes near free TE: {len(nongenes_near_free_TE_unique)} / {len(nonretrogenes)} ({len(nongenes_near_free_TE_unique)/len(nonretrogenes):.2%})")


	table = [
    [len(rgenes_near_free_TE_unique), len(retrogenes)-len(rgenes_near_free_TE_unique)],     # Retrogenes
    [len(nongenes_near_free_TE_unique), len(nonretrogenes)-len(nongenes_near_free_TE_unique)]    # Non-retrogenes
	]

	chi2, p, dof, expected = chi2_contingency(table)

	rg_te_proximity_vec.append(len(rgenes_near_free_TE_unique)/len(retrogenes))
	nrg_te_proximity_vec.append(len(nongenes_near_free_TE_unique)/len(nonretrogenes))
	P_val_list.append(p)


print(rg_te_proximity_vec)
print(nrg_te_proximity_vec)
print(P_val_list)

# === Plot ===
prox_list = [x - 1 for x in prox_list]

# Compute -log10(p-values)
log_pvals = [-np.log10(p) if p > 0 else 300 for p in P_val_list]

def get_significance_star(p):
    """Return number of stars based on p-value thresholds."""
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''

def plot_proximity_with_significance(prox_list, rg_te_proximity_vec, nrg_te_proximity_vec, p_val_list,
                                     title="Proximity to retrotransposons", ylabel="Fraction of genes near retrotransposons"):
    """
    Plot proximity curve for retrogenes and non-retrogenes, with significance stars and legend.

    Parameters:
    - prox_list: list of distance thresholds (x-axis)
    - rg_te_proximity_vec: list of retrogene fractions
    - nrg_te_proximity_vec: list of non-retrogene fractions
    - p_val_list: list of p-values (for chi-square tests at each proximity)
    """
    plt.figure(figsize=(8, 5))

    # Plot lines
    plt.plot(prox_list, rg_te_proximity_vec, marker='o', label='Retrogenes', linewidth=2)
    plt.plot(prox_list, nrg_te_proximity_vec, marker='s', label='Non-retrogenes', linewidth=2)

    # Annotate significance stars
    for x, y1, y2, p in zip(prox_list, rg_te_proximity_vec, nrg_te_proximity_vec, p_val_list):
        star = get_significance_star(p)
        if star:
            ymax = max(y1, y2)
            plt.text(x, ymax + 0.02, star, ha='center', va='bottom', fontsize=12, color='black')

    # Axis and labels
    plt.xlabel("Proximity to retrotransposons (bp)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(prox_list, rotation=45)
    plt.grid(False)

    # Legend with significance explanation
    base_handles, base_labels = plt.gca().get_legend_handles_labels()
    star_legend = [
        Line2D([0], [0], color='none', marker='', label='*   p < 0.05'),
        Line2D([0], [0], color='none', marker='', label='**  p < 0.01'),
        Line2D([0], [0], color='none', marker='', label='*** p < 0.001')
    ]
    plt.legend(base_handles + star_legend, base_labels + [h.get_label() for h in star_legend], loc='upper left')

    plt.tight_layout()
    plt.show()


# Example usage (assuming you’ve already calculated these lists):
plot_proximity_with_significance(
    prox_list,
    rg_te_proximity_vec,
    nrg_te_proximity_vec,
    P_val_list
)

