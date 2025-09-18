#!/bin/bash
#SBATCH --job-name=atac_multicov
#SBATCH --time=5:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4

# Load bedtools if needed
module load BEDTools

cd /home/zl436/palmer_scratch/bam/ATAC 
bedtools multicov -bams *.bam -bed consensus_peaks_min3_merged.bed > peak_counts_consensus.txt

