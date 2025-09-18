#!/bin/bash
#SBATCH --job-name=extract_variants     # Job name
#SBATCH --output=extract_variants_%j.log  # Standard output and error log
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=4G                        # Total memory required (adjust as needed)
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec


module load VCFtools

# Submit the job (example of VCF, chromosome, start, end, and output file)
VCF_FILE="/home/zl436/ycga_work/oyster_genome/vcf/20cv_NC_035784.1.vcf"
CHROM="NC_035784.1"
START="58754900"
END="58755600"
OUTPUT_FILE="/home/zl436/ycga_work/oyster_genome/vcf/20cv_${CHROM}_${START}_${END}.vcf"

vcftools --vcf ${VCF_FILE} --chr ${CHROM} --from-bp ${START} --to-bp ${END} --recode --out ${OUTPUT_FILE}

