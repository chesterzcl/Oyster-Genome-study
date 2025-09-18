#!/bin/bash
#SBATCH --job-name=vcf_maf_filter
#SBATCH --time=03:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=2

# Load modules if needed
module load BCFtools tabix

# Input and output paths
VCF_IN="/home/zl436/palmer_scratch/vcf/oyster/snp/2nd_draft_new/20cv_df3_hoh_biSNP_filtered.vcf"
VCF_OUT="/home/zl436/palmer_scratch/vcf/oyster/snp/2nd_draft_new/filtered_snps_maf05.vcf.gz"

bcftools +fill-tags "$VCF_IN" -- -t MAF | \
bcftools view -i 'MAF>=0.05' -Oz -o "$VCF_OUT"

# Step 3: Index output VCF
tabix -p vcf "$VCF_OUT"


