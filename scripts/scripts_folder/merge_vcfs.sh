#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL

module load Java
module load GATK

dir=/home/zl436/palmer_scratch/vcf/oyster/2nd_draft

gatk MergeVcfs -I ${dir}/20cv_df2_hoh_HiC_scaffold_1_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_2_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_3_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_4_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_5_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_6_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_7_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_8_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_9_biSNP_filtered.vcf -I ${dir}/20cv_df2_hoh_HiC_scaffold_10_biSNP_filtered.vcf -O ${dir}/20cv_df2_hoh_biSNP_filtered.vcf
