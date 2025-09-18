#!/bin/bash
#SBATCH --job-name=cluster_cnv
#SBATCH --time=02:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4

module load BCFtools
module load BEDTools

# === Directories ===
VCF_DIR="filtered_vcfs"
BED_DIR="clustered_beds"
mkdir -p $BED_DIR

RAW_BED="$BED_DIR/all_samples_raw.bed"
CLUSTERED="$BED_DIR/all_samples_clustered.bed"

> "$RAW_BED"

# === Extract CNV calls into BED format with CNV type and sample name ===
for f in $VCF_DIR/*.vcf.gz; do
    sample=$(basename "$f" .filtered.vcf.gz)

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ALT\n' "$f" \
      | sed 's/<//; s/>//' | awk -v OFS='\t' -v s="$sample" '{print $1, $2, $3, $4, s}' \
      >> "$RAW_BED"

done

# === Cluster separately by CNV type ===
> "$CLUSTERED"

for svtype in DEL DUP; do
    awk -v type="$svtype" '$4 == type' "$RAW_BED" | \
    bedtools sort -i - | \
    bedtools cluster -d 1000 -i - | \
    awk -v OFS='\t' -v t="$svtype" '{$4 = t; print}' >> "$CLUSTERED"
done

module load miniconda

conda activate seq_env

python make_cnv_matrix_from_clustered.py


