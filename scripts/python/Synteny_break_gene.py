import pandas as pd
import pybedtools
import os 

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final")
# ================================
# CONFIG
# ================================
FLANK_BED = "breakpoints_cf_inside_gap_windows.bed"
GENE_BED = "genes_function.bed"
OUTPUT_OVERLAP = "breakend_flank_gene_overlaps.csv"

# ================================
# 1. Load flank windows and gene BED as BedTool
# ================================
flank_bt = pybedtools.BedTool(FLANK_BED)
gene_bt = pybedtools.BedTool(GENE_BED)

# ================================
# 2. Intersect
# ================================
# flank fields: 0-chrom, 1-start, 2-end, 3-name, 4-score, 5-event_type
# gene fields: 6-chrom, 7-start, 8-end, 9-gene_id, 10-gene_name, 11-strand
overlaps = flank_bt.intersect(gene_bt, wa=True, wb=True)

# ================================
# 3. Parse results
# ================================
overlap_records = []
for f in overlaps:
    fields = f.fields
    overlap_records.append({
        "flank_chrom": fields[0],
        "flank_start": int(fields[1]),
        "flank_end": int(fields[2]),
        "flank_name": fields[3],
        "event_type": fields[5],
        "gene_chrom": fields[6],
        "gene_start": int(fields[7]),
        "gene_end": int(fields[8]),
        "gene_id": fields[9],
        "gene_name": fields[10] if fields[10] else None,
        "strand": fields[11]
    })

overlap_df = pd.DataFrame(overlap_records)

print(f"Found {len(overlap_df)} overlapping gene-flank pairs")

# ================================
# 4. Add flank direction (parsed from name)
# ================================
overlap_df["flank_side"] = overlap_df["flank_name"].apply(
    lambda x: "downstream" if "downstream" in x else "upstream"
)

# ================================
# 5. Save CSV
# ================================
overlap_df.to_csv(OUTPUT_OVERLAP, index=False)
print(f"*** Saved {OUTPUT_OVERLAP} ***")


# ================================
# 6. Generate foreground and background gene NAME lists
# ================================
print("\n*** STEP 6: Generate foreground and background gene NAME lists for GO analysis ***")

# 6.1 Load genes_function.bed for background
gene_cols = ["chrom","start","end","gene_id","gene_name","strand"]
gene_df = pd.read_csv(GENE_BED, sep="\t", header=None, names=gene_cols)

# Filter to those with non-empty gene_name
background_genes = gene_df["gene_name"].dropna().str.strip()
background_genes = background_genes[background_genes != ""].unique()
print(f"Extracted {len(background_genes)} unique background gene NAMES with annotation")

with open("background_gene_name_list.txt", "w") as f:
    for gene in sorted(background_genes):
        f.write(f"{gene}\n")
print("*** Saved background_gene_name_list.txt ***")

# 6.2 Extract foreground gene names from overlaps
if "gene_name" in overlap_df.columns:
    foreground_genes = overlap_df["gene_name"].dropna().str.strip()
    foreground_genes = foreground_genes[foreground_genes != ""].unique()
    print(f"Extracted {len(foreground_genes)} unique foreground gene NAMES with annotation")

    with open("foreground_gene_name_list.txt", "w") as f:
        for gene in sorted(foreground_genes):
            f.write(f"{gene}\n")
    print("*** Saved foreground_gene_name_list.txt ***")
else:
    print("WARNING: No gene_name column found in overlap_df; skipping foreground list generation.")

# ================================
# 7. Perform GO enrichment with gprofiler
# ================================
from gprofiler import GProfiler

print("\n*** STEP 7: Run GO enrichment with gprofiler ***")

# Load foreground and background gene name lists from files
with open("foreground_gene_name_list.txt") as f:
    foreground_genes = [line.strip() for line in f if line.strip()]

with open("background_gene_name_list.txt") as f:
    background_genes = [line.strip() for line in f if line.strip()]

print(f"Foreground set size: {len(foreground_genes)}")
print(f"Background set size: {len(background_genes)}")

# Initialize g:Profiler
gp = GProfiler(return_dataframe=True)

# Run enrichment
# Organism: e.g. mollusk is likely not in the standard list, so let's use "mus_musculus" or "drosophila_melanogaster" as placeholder
# Replace with appropriate organism if you have ortholog mapping!
result_df = gp.profile(
    organism='hsapiens',
    query=foreground_genes,
    user_threshold=0.05,
    significance_threshold_method="fdr",
    # background=background_genes,
    sources=["GO:BP"],
    no_evidences=True
)

# Save result
result_df.to_csv("gprofiler_enrichment_results.csv", index=False)
print("*** Saved gprofiler_enrichment_results.csv ***")

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ================================
# 1. Filter significant GO:BP results
# ================================
sig_results = result_df[result_df['p_value'] < 0.05].copy()

if sig_results.empty:
    print("*** No significant terms at FDR < 0.05 ***")
else:
    print(f"*** Found {len(sig_results)} significant terms ***")

    # ================================
    # 2. Compute enrichment ratio
    # ================================
    # enrichment_ratio = (query_gene_count / query_size) / (term_gene_count / term_size)
    sig_results['enrichment_ratio'] = (
        (sig_results['intersection_size'] / sig_results['query_size']) /
        (sig_results['term_size'] / sig_results['effective_domain_size'])
    )

    # ================================
    # 3. Plot top terms
    # ================================
    top_n = 10
    plot_df = sig_results.nsmallest(top_n, 'p_value').copy()
    plot_df = plot_df.sort_values('p_value', ascending=True)

    plt.figure(figsize=(8, 0.5 * len(plot_df) + 4))
    barplot = sns.barplot(
        data=plot_df,
        x='enrichment_ratio',
        y='name',
        hue=-np.log10(plot_df['p_value']),
        palette='viridis'
    )
    barplot.set_xlabel('Enrichment Ratio',fontsize=13)
    barplot.set_ylabel('GO Biological Process Term',fontsize=13)
    plt.legend(title='-log10(p-value)')
    # plt.title('Top Enriched GO:BP Terms')
    plt.tight_layout()
    plt.savefig("gprofiler_enrichment_top_terms.png")
    plt.close()

    print("*** Saved gprofiler_enrichment_top_terms.png ***")

