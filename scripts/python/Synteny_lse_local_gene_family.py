import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

def plot_gene_cluster_with_TEs(gff_file, repeat_bed, chrom, start, end, output_file):
    # --- Load gene annotation ---
    gff_cols = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff = pd.read_csv(gff_file, sep='\t', comment='#', names=gff_cols)
    gff = gff[(gff['seqid'] == chrom) & (gff['end'] >= start) & (gff['start'] <= end)]

    # --- Extract gene names ---
    def extract_gene_name(attr):
        for part in attr.split(";"):
            if part.startswith("Name="):
                return part.split("=")[1]
        return None

    gff['gene_name'] = gff['attributes'].apply(extract_gene_name)
    gff['gene_id'] = gff['attributes'].str.extract(r'ID=([^;]+)')

    # Step 1: Build mapping from gene_id to longest mRNA ID
    gene_to_longest_mrna = {}

    for gene_id in gff[gff['type'] == 'gene']['gene_id'].unique():
        mRNAs = gff[(gff['type'] == 'mRNA') & (gff['attributes'].str.contains(f'Parent={gene_id}'))].copy()
        if not mRNAs.empty:
            mRNAs.loc[:, 'length'] = mRNAs['end'] - mRNAs['start']
            longest_row = mRNAs.loc[mRNAs['length'].idxmax()]
            mrna_id = longest_row['attributes'].split(";")[0].split("=")[1]
            gene_to_longest_mrna[gene_id] = mrna_id

    # Step 2: Filter gff to keep only features belonging to longest mRNA
    def belongs_to_longest_mrna(row):
        for gene_id, mrna_id in gene_to_longest_mrna.items():
            if f"Parent={mrna_id}" in row['attributes']:
                return True
        return row['type'] in ['gene', 'mRNA']  # Keep top-level too

    gff = gff[gff.apply(belongs_to_longest_mrna, axis=1)]

    # --- Track assignment ---

    gene_tracks = {}
    for gene_id, group in gff[gff['type'] == 'gene'].groupby('gene_id'):
        name = group['gene_name'].iloc[0]
        if name:
            gene_tracks[gene_id] = 1  # Named genes on track 1
        else:
            gene_tracks[gene_id] = 2  # All unnamed genes share track 2

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(12, 3))
    color_map = {'+': '#4daf4a', '-': '#e41a1c'}

    # Build gene_id ➝ mRNA_id map
    gene2mrna = {}
    for _, row in gff[gff['type'] == 'mRNA'].iterrows():
        gene_id = row['attributes'].split("Parent=")[-1].split(";")[0]
        mrna_id = row['attributes'].split("ID=")[-1].split(";")[0]
        gene2mrna.setdefault(gene_id, []).append(mrna_id)

    # Loop over gene_id from type == gene
    for gene_id in gff[gff['type'] == 'gene']['gene_id'].unique():
        track = gene_tracks.get(gene_id, 1)
        color = color_map.get(gff[gff['gene_id'] == gene_id]['strand'].iloc[0], 'gray')
        
        # For each transcript
        for mrna_id in gene2mrna.get(gene_id, []):
            gene_df = gff[gff['attributes'].str.contains(f'Parent={mrna_id}')]
            # Plot exons
            for _, row in gene_df[gene_df['type'] == 'exon'].iterrows():

                ax.add_patch(patches.Rectangle(
                    (row['start'], track - 0.2),
                    row['end'] - row['start'],
                    0.4,
                    facecolor=color,
                    edgecolor=color,
                    linewidth=0.3
                ))

            # Determine the min/max positions of exons
            exons = gene_df[gene_df['type'] == 'exon']
            if not exons.empty:
                exon_start = exons['start'].min()
                exon_end = exons['end'].max()
                ax.hlines(y=track, xmin=exon_start, xmax=exon_end,
                          color="black", linewidth=0.5, linestyle='-')
        
        # Add label
        gene_mid = gff[gff['gene_id'] == gene_id][['start', 'end']].agg(['min', 'max']).mean().mean()
        label = gene_id
        ax.text(gene_mid, track + 0.3, label, ha='center', fontsize=7)

    # --- Plot TEs ---
    repeat_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'class']
    try:
        rep = pd.read_csv(repeat_bed, sep='\t', header=None)
        rep.columns = repeat_cols[:rep.shape[1]]  # allow for missing columns
    except:
        rep = pd.read_csv(repeat_bed, sep='\t', header=None, names=repeat_cols)

    rep = rep[(rep['chrom'] == chrom) & (rep['start'] >= start) & (rep['end'] <= end)]
    rep['class'] = rep['class'].fillna('Unknown')

    # --- TE Class-specific Tracks ---
    te_colors = {'LINE': 'blue', 'SINE': 'orange', 'DNA': 'purple', 'LTR': 'green'}
    te_tracks = {'LINE': -0.2, 'SINE': -0.6, 'DNA': -1.0, 'LTR': -1.4}  # Y positions per class

    for _, row in rep.iterrows():
        matched_class = next((k for k in te_colors if k in str(row['class'])), 'Unknown')
        color = te_colors.get(matched_class, 'gray')
        y_track = te_tracks.get(matched_class, -1.8)  # put unknowns at bottom

        ax.add_patch(patches.Rectangle(
            (row['start'], y_track), 
            row['end'] - row['start'], 
            0.3,
            facecolor=color, alpha=0.6
        ))

    # Extend Y-axis limits to include TE tracks
    ax.set_ylim(min(te_tracks.values()), 3)  # Only 2 gene tracks now

    # --- Formatting ---
    # Format x-axis to show Mb
    xticks = ax.get_xticks()
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{x/1e6:.2f}" for x in xticks],fontsize=9)
    # ax.set_xlabel(f"Genomic position on {chrom}(mb)")
    ax.set_ylabel("Genomic Elements")
    # ax.set_title(f"Gene architecture and TEs from {start:,} to {end:,}")
    # --- Y-axis Labels ---
    yticks = [1, 2]  # Track 1: Named genes, Track 2: Unnamed genes
    ylabels = ["FRRS1 Genes", "Other Genes"]

    # Add TE track labels (match te_tracks)
    te_tracks = {'LINE': -0.2, 'SINE': -0.6, 'DNA': -1.0, 'LTR': -1.4}

    # Add TE tracks centered slightly higher
    yticks += [y + 0.15 for y in te_tracks.values()]
    ylabels += list(te_tracks.keys())
    
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels,fontsize=9)
    # Only keep x and y axes, remove top and right
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xlim(start, end)  # Enforce strict x-limits
    ax.set_autoscale_on(False)  # Disable automatic axis expansion
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

# === Load GFF ===
gff_file = "braker_with_name_and_description.gff3"  # Update to your actual path

buffer=10000
chr="HiC_scaffold_1"
start=13283141
end=13446823

plot_gene_cluster_with_TEs(
    gff_file=gff_file,
    repeat_bed="all_repeats.bed",
    chrom=chr,
    start=start-buffer,
    end=end+buffer,
    output_file="FRRS1_cluster_with_TEs.png"
)






