#!/bin/bash
#SBATCH --job-name=merge_vcfs
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G

# Load GATK module (adjust if needed)
module load GATK

# Output file name
OUTPUT="20cv_df3_hoh.vcf"

# Run GATK GatherVcfs
gatk GatherVcfs \
  -I 20cv_df3_hoh_HiC_scaffold_1.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_2.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_3.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_4.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_5.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_6.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_7.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_8.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_9.vcf \
  -I 20cv_df3_hoh_HiC_scaffold_10.vcf \
  -O "$OUTPUT"

# Index the resulting VCF
gatk IndexFeatureFile -I "$OUTPUT"
