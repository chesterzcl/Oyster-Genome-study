import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches


fai_file ="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/primary_dedup_chr.fa.fai"
gene_file="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/genes_function.bed"
TE_dna_file="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/repeat_analysis/dna_transposons.bed"
TE_retro_file="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/repeat_analysis/retroelements.bed"
stlt_file="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/satellite_analysis_satellite_filtered.bed"
out_dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/"

# === Load chromosome sizes ===
fai = pd.read_csv(fai_file,sep="\t",header=None,usecols=[0,1],names=["chrom", "length"])
fai = fai[["chrom", "length"]]

gene=pd.read_csv(gene_file,sep='\t',header=None,names=["chrom","start","end","id","gene_name","subclass"])
gene=gene[gene["subclass"]!="unknown"]
gene["type"]="Protein-coding gene"
print(gene)

TE_dna=pd.read_csv(TE_dna_file,sep='\t',header=None,usecols=[0,1,2,3],names=["chrom","start","end","id"])
TE_retro=pd.read_csv(TE_retro_file,sep='\t',header=None,usecols=[0,1,2,3],names=["chrom","start","end","id"])
TE_dna["type"]="DNA transposon"
TE_retro["type"]="Retrotransposon"
TE=pd.concat([TE_dna,TE_retro],axis=0,ignore_index=True)
print(TE)

stlt=pd.read_csv(stlt_file,sep='\t',header=None,usecols=[0,1,2,3],names=["chrom","start","end","id"])
stlt["type"]="Satellite"
print(stlt)

combined = pd.concat([
    gene[['chrom','start','end','type']],
    TE_dna[['chrom','start','end','type']],
    TE_retro[['chrom','start','end','type']],
    stlt[['chrom','start','end','type']]
],ignore_index=True)


combined['chrom_num']=combined['chrom'].str.extract(r'_(\d+)$').astype(int)
combined = combined.sort_values(['chrom_num','start']).drop(columns='chrom_num').reset_index(drop=True)
combined.to_csv(f"{out_dir}all_features.bed",sep="\t",index=False,header=False)

print(combined)

fai = fai.copy()
fai['chrom_num'] = fai['chrom'].str.extract(r'_(\d+)$').astype(int)
fai = fai.sort_values('length',ascending=False).reset_index(drop=True)

# Color map for feature types
colors = {
    "Protein-coding gene": "#4daf4a",
    "DNA transposon":      "#e41a1c",
    "Retrotransposon":     "#377eb8",
    "Satellite":           "#984ea3"
}

# Plot setup
fig,ax=plt.subplots(figsize=(12,len(fai)*0.4))  
chrom_height=0.5
feature_height=chrom_height*0.95  # feature bars slightly thinner than chromosome
overlap = (chrom_height-feature_height)/2  # vertical offset for centering feature bars

for s in ['top','left','right']:
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_linewidth(1)

# Adjust margins to prevent overlap
plt.subplots_adjust(left=0.95,right=1,top=0.5,bottom=0.15)

for i, row in fai.iterrows():
    chrom=row['chrom']
    length_mb=row['length']/1e6
    chrom_idx=row['chrom_num']
    y=i

    # Draw background chromosome bar with thinner black boundary
    ax.add_patch(
        patches.Rectangle(
            (0, y), length_mb, chrom_height,
            facecolor='lightgray',edgecolor='black',linewidth=0.4
        )
    )

    # Overlay each feature segment with no border, centered vertically
    sub = combined[combined['chrom'] == chrom]
    for _, feat in sub.iterrows():
        start_mb=feat['start']/1e6
        end_mb=feat['end']/1e6
        width=end_mb-start_mb
        ax.add_patch(
            patches.Rectangle(
                (start_mb,y+overlap),width,feature_height,
                facecolor=colors.get(feat['type'],'gray'),
                edgecolor='none',
                linewidth=0.1
            )
        )

    # Label with chrX format
    label=f"chr{chrom_idx}"
    ax.text(-0.05,y+chrom_height/2,label,va='center',ha='right',fontsize=7)

# Styling axes
ax.set_xlim(-2,fai['length'].max()/1e6+2)
ax.set_ylim(-0.5,len(fai)-0.5)
ax.set_yticks([])
ax.set_xlabel("Position (Mb)")
ax.set_title("Chromosome-level Distribution of Genome Features")

# Legend with smaller markers and better placement
handles=[patches.Patch(color=c,label=k) for k,c in colors.items()]
ax.legend(handles=handles,loc='upper right',bbox_to_anchor=(0.98,1),borderaxespad=0.5,frameon=False,fontsize=7)

plt.tight_layout()
plt.savefig(f"{out_dir}chromosome_feature_distribution.png",dpi=300)






