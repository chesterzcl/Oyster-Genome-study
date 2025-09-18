#!/bin/bash
#SBATCH --job-name=plink_rel_pca
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# === Load PLINK 1.9b ===
module load PLINK/1.9b_6.21-x86_64

# === Config ===
VCF=/home/zl436/palmer_scratch/vcf/oyster/snp/2nd_draft_new/20cv_df3_hoh_biSNP_filtered_maf05.vcf.gz
OUT_DIR=/home/zl436/palmer_scratch/script/ld_analysis
OUT_PREFIX=cv20_maf05

mkdir -p $OUT_DIR

# === Step 1: Convert VCF to PLINK binary format ===
#plink --vcf $VCF \
#      --make-bed \
#      --double-id \
#      --allow-extra-chr \
#      --out $OUT_DIR/${OUT_PREFIX}

# === Step 2: Calculate pairwise relatedness (PI_HAT) ===
plink --bfile $OUT_DIR/${OUT_PREFIX} \
      --genome full \
      --allow-extra-chr \
      --out $OUT_DIR/${OUT_PREFIX}_relatedness \
      --threads 4

# === Step 3: Principal Component Analysis (PCA) ===
plink --bfile $OUT_DIR/${OUT_PREFIX} \
      --allow-extra-chr \
      --pca 10 \
      --out $OUT_DIR/${OUT_PREFIX}_pca \
      --threads 4

echo "✅ Relatedness and PCA completed for ${OUT_PREFIX}"


