#!/bin/bash
#SBATCH --job-name=ld_decay
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8

module load PLINK/1.9b_6.21-x86_64

# === Config ===
VCF=/home/zl436/palmer_scratch/vcf/oyster/snp/2nd_draft_new/20cv_df3_hoh_biSNP_filtered_maf05.vcf.gz
OUT_DIR=/home/zl436/palmer_scratch/script/ld_analysis
mkdir -p $OUT_DIR
THREADS=8

mkdir -p $OUT_DIR

# Convert VCF to PLINK format (treat HiC_scaffold_X as chromosomes)
plink --vcf $VCF \
      --double-id \
      --allow-extra-chr \
      --make-bed \
      --out $OUT_DIR/oyster_all

# Thin SNPs to reduce computation (optional but recommended for LD decay)
# This keeps ~8% of SNPs randomly (~1 per 200bp if density is high)
plink --bfile $OUT_DIR/oyster_all \
      --thin 0.1 \
      --allow-extra-chr \
      --make-bed \
      --out $OUT_DIR/oyster_thinned

# Compute pairwise LD (r²) within 500kb windows across all scaffolds
plink --bfile $OUT_DIR/oyster_thinned \
      --r2 \
      --ld-window-kb 500 \
      --ld-window 500 \
      --ld-window-r2 0.2 \
      --allow-extra-chr \
      --threads 8 \
      --out $OUT_DIR/ld_500kb_r2

# Completion message
echo "LD calculation done. Output: $OUT_DIR/ld_500kb_r2.ld"

