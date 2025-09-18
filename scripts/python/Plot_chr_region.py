import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly")
# --- Inputs ---
gff3="braker.gff3"



region = ("HiC_scaffold_1",13258944,13446823)

# --- Load gene models from GFF3 ---
cols = ["chrom","source","feature","start","end","score","strand","phase","attrs"]
gff = pd.read_csv(
    gff3,
    sep="\t",
    comment="#",
    header=None,
    names=cols,
    dtype={"attrs":str}
)


# Keep only features in region
chrom, win_s, win_e = region
g = gff[(gff.chrom == chrom) & (gff.end >= win_s) & (gff.start <= win_e)].copy()




# Extract gene and transcript IDs
g["gene_id"] = g['attrs'].str.extract(r"Parent=([^;]+)")
g.loc[g.feature=="gene", "gene_id"] = g.loc[g.feature=="gene", "attrs"].str.extract(r"ID=([^;]+)")

# Shift to window-relative coords
g["s2"] = g.start.clip(lower=win_s) - win_s
g["e2"] = g.end.clip(upper=win_e) - win_s


# Prepare plotting tracks per transcript
tx_ids = g[g.feature=="mRNA"]["attrs"].str.extract(r"ID=([^;]+)").iloc[:,0].unique()
n = len(tx_ids)+3
track_height = 0.2
track_gap    = 0.2

combined = pd.read_csv(f"all_features.bed", sep="\t", names=["chrom","start","end","type"])
sub_te = combined[(combined['chrom']==chrom) &
                  (combined['end'] >= win_s) & (combined['start'] <= win_e)].copy()
sub_te['s2'] = sub_te['start'].clip(lower=win_s) - win_s
sub_te['e2'] = sub_te['end'].clip(upper=win_e)   - win_s
# --- Plot settings ---
track_types = [
    ('Satellite', 'Satellite', 0.2),
    ('Retrotransposon', 'Retrotransposon', 0.2),
    ('DNA transposon', 'DNA transposon', 0.2)
]

# then transcript tracks
for tx in tx_ids:
    track_types.append(('Gene', tx, 0.3))

# colors
domain_colors = {
    'DNA transposon':      '#e41a1c',
    'Retrotransposon':     '#377eb8',
    'Satellite':           '#984ea3'
}

fig,ax=plt.subplots(figsize=(8, n*(track_height+track_gap)+1))

# For each transcript, plot introns (lines) and exons (blocks), label gene
y=n

for i, tx in enumerate(tx_ids):
    sub = g[g['attrs'].str.contains(f"Parent={tx}")]
    # invert so first on top
    y=y-1
    # Label with gene name (use gene_id of transcript)
    gene = g[g['feature']=="mRNA"][g['attrs'].str.contains(f"ID={tx}")]["gene_id"].iloc[0]
    ax.text(- (win_e-win_s)*0.01, y*(track_height+track_gap)+track_height/2,
            gene, va="center", ha="right", fontsize=7)

    # Draw intron as a horizontal line between transcript start/end
    tx_start = sub["s2"].min()
    tx_end   = sub["e2"].max()
    ax.hlines(y*(track_height+track_gap)+track_height/2,
              tx_start, tx_end, color="black", linewidth=0.5)

    # Draw exons
    for _, exon in sub[sub.feature=="exon"].iterrows():
        ax.add_patch(patches.Rectangle(
            (exon.s2, y*(track_height+track_gap)), exon.e2-exon.s2,track_height,
            facecolor="#4daf4a", edgecolor="black", linewidth=0.1
        ))
    

# Draw repeats

for repeat_type in domain_colors.keys():
    y=y-1
    ax.text(- (win_e-win_s)*0.01, y*(track_height+track_gap)+track_height/2,
                repeat_type, va="center", ha="right", fontsize=7)   
    for _,repeat in sub_te[sub_te["type"]==repeat_type].iterrows():
        ax.add_patch(patches.Rectangle(
                (repeat.s2,y*(track_height+track_gap)),repeat.e2-repeat.s2,track_height,
                facecolor=domain_colors[repeat_type],edgecolor="black",linewidth=0.1
                ))

#Vertical lines
plt.axvline(x=50101934-win_s, color='b', linestyle='--',linewidth=0.5)
plt.axvline(x=50102809-win_s, color='b', linestyle='--',linewidth=0.5)

# Styling
ax.set_xlim(0, win_e-win_s)
# Convert x-axis to absolute coordinates
xt = ax.get_xticks()
# compute absolute positions and set labels
abs_xt = (xt + win_s).astype(int)
ax.set_xticks(xt)
ax.set_xticklabels([f"{x:,}" for x in abs_xt], rotation=0, fontsize=7)
ax.set_xlabel(f"{chrom}:{win_s:,}-{win_e:,}")
ax.set_ylim(0,n/2)
ax.set_yticks([])
for spine in ['top','right','left']:
    ax.spines[spine].set_visible(False)
ax.spines['bottom'].set_linewidth(0.8)
ax.set_xlabel(f"Position on Chr{region[0].split('_')[-1]}")


plt.tight_layout()
plt.savefig(f"genome_track_Chr{region[0].split('_')[-1]}_{region[1]}_{region[2]}",dpi=300)


