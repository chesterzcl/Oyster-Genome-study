#!/bin/bash
#SBATCH --job-name=snp_het_cluster
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

module load BCFtools
module load BEDTools

# === Paths ===
VCF=/home/zl436/palmer_scratch/vcf/oyster/snp/2nd_draft_new/20cv_df3_hoh_biSNP_filtered.vcf
BED_DIR=/home/zl436/palmer_scratch/bed
OUT_DIR=/home/zl436/palmer_scratch/ATAC_analysis
mkdir -p $OUT_DIR

# === Step 1: Calculate heterozygosity for each SNP ===
echo "[INFO] Calculating SNP heterozygosity..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $VCF \
| awk '
BEGIN { OFS="\t" }
{
  het=0; total=0
  for(i=5; i<=NF; i++) {
    if($i=="0/1" || $i=="1/0" || $i=="0|1" || $i=="1|0") het++  
    if($i!="." && $i!="./.") total++
  }
  if (total > 16) {
    print $1, $2-1, $2, het/total
  }
}' > $OUT_DIR/site_het.bed

# === Step 2: Intersections ===
echo "[INFO] Intersecting with genomic annotations..."

# Genes
bedtools intersect -a $OUT_DIR/site_het.bed \
                   -b $BED_DIR/genes.bed \
                   -wa > $OUT_DIR/het_in_genes.bed

# 2kb Promoters
bedtools intersect -a $OUT_DIR/site_het.bed \
                   -b $BED_DIR/promoters_2kb.bed \
                   -wa > $OUT_DIR/het_in_promoters2kb.bed

# High Variance ATAC Peaks
bedtools intersect -a $OUT_DIR/site_het.bed \
                   -b $BED_DIR/ATAC_high_variance_peaks.bed \
                   -wa > $OUT_DIR/het_in_ATAC_high.bed

# Low Variance ATAC Peaks
bedtools intersect -a $OUT_DIR/site_het.bed \
                   -b $BED_DIR/ATAC_low_variance_peaks.bed \
                   -wa > $OUT_DIR/het_in_ATAC_low.bed

# High Variance Peaks Overlapping Promoters
bedtools intersect -a $OUT_DIR/site_het.bed \
                   -b $BED_DIR/ATAC_high_variance_peaks_promoter_annotated.bed \
                   -wa > $OUT_DIR/het_in_ATAC_high_promoter.bed

# Low Variance Peaks Overlapping Promoters
bedtools intersect -a $OUT_DIR/site_het.bed \
                   -b $BED_DIR/ATAC_low_variance_peaks_promoter_annotated.bed \
                   -wa > $OUT_DIR/het_in_ATAC_low_promoter.bed


# === Completion Message ===
echo "[DONE] Output files written to $OUT_DIR:"
ls -lh $OUT_DIR/*.bed


