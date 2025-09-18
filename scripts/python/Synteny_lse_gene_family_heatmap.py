import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/MG_CV_final/wgd_ksd_FRRS1")

def plot_ka_ks_heatmap_with_two_colorbars(wgd_file, output_file=None):
    # Load data
    df = pd.read_csv(wgd_file, sep='\t')
    df['gene1_short'] = df['gene1'].str.extract(r'^(g\d+)')
    df['gene2_short'] = df['gene2'].str.extract(r'^(g\d+)')

    gene_ids = sorted(set(df['gene1_short']).union(set(df['gene2_short'])))
    n = len(gene_ids)
    gene_idx = {gene: i for i, gene in enumerate(gene_ids)}

    # Initialize matrices
    ks_mat = np.full((n, n), np.nan)
    omega_mat = np.full((n, n), np.nan)

    for _, row in df.iterrows():
        i, j = gene_idx[row['gene1_short']], gene_idx[row['gene2_short']]
        if i != j:
            ks_mat[i, j] = ks_mat[j, i] = row['dS']
            omega_mat[i, j] = omega_mat[j, i] = row['dN/dS']

    # Create figure and gridspec
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[20, 1], width_ratios=[20, 1], hspace=0.15, wspace=0.15)

    # Main heatmap axes
    ax = fig.add_subplot(gs[0, 0])
    mask_upper = np.triu(np.ones((n, n), dtype=bool))
    mask_lower = np.tril(np.ones((n, n), dtype=bool))

    sns.heatmap(ks_mat, mask=mask_upper, cmap="Blues", ax=ax,
                cbar=False, annot=True, fmt=".2f", linewidths=0.5, linecolor='white')
    sns.heatmap(omega_mat, mask=mask_lower, cmap="Reds", ax=ax,
                cbar=False, annot=True, fmt=".2f", linewidths=0.5, linecolor='white')

    ax.set_aspect('equal')
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_xticklabels(gene_ids)
    ax.set_yticklabels(gene_ids)
    # ax.set_title("Pairwise Ka/Ks Heatmap")

    # ω (Ka/Ks) colorbar (vertical, shrunk to 60%)
    sm1 = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(vmin=np.nanmin(omega_mat), vmax=np.nanmax(omega_mat)))
    sm1.set_array([])
    cbar1 = plt.colorbar(sm1, ax=ax, shrink=0.3, location='right')  # auto-position right of heatmap
    cbar1.set_label('ω (Ka/Ks)', rotation=90)

    # Ks colorbar (horizontal, shrunk to 70%)
    sm2 = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin=np.nanmin(ks_mat), vmax=np.nanmax(ks_mat)))
    sm2.set_array([])
    cbar2 = plt.colorbar(sm2, ax=ax, orientation='horizontal', shrink=0.4, pad=0.15)
    cbar2.set_label('Ks')

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    else:
        plt.tight_layout()
        plt.show()


# Example usage
plot_ka_ks_heatmap_with_two_colorbars("FRRS1.fa.tsv.ks.tsv", output_file="ka_ks_heatmap.png")

