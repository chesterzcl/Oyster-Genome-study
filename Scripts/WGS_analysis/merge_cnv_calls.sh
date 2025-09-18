#!/bin/bash

#SBATCH --job-name=filter_merge_cnv
#SBATCH --time=02:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=2

module load BCFtools

# Create output directory
OUTPUT_DIR="filtered_vcfs"
mkdir -p $OUTPUT_DIR

echo "Filtering CNV VCFs..."
for f in genotyped-segments-*.vcf.gz; do
    sample=$(basename "$f" .vcf.gz)
    outvcf="${OUTPUT_DIR}/${sample}.filtered.vcf.gz"

    # Filter CNV calls
    bcftools view \
        -i 'GT!="0/0" && CN!=2 && QS>=30 && QA>=20 && QSE>=20 && QSS>=20' \
        "$f" -Oz -o "$outvcf"

    bcftools index "$outvcf"
done

echo "Merging filtered CNV VCFs..."
FILTERED_VCF_LIST=$(ls $OUTPUT_DIR/*.filtered.vcf.gz)

# Merge all filtered VCFs into one multisample VCF
bcftools merge $FILTERED_VCF_LIST -Oz -o merged_filtered_cnv.vcf.gz

# Index final merged VCF
bcftools index merged_filtered_cnv.vcf.gz

echo "Done."


