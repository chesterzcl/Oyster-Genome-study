#!/bin/bash
#SBATCH --job-name=consensus_peaks
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4


###########################################################
# Script: extract_fragment_lengths.sh
# Purpose: Extract fragment length (TLEN) distributions 
#          from paired-end BAM files using samtools
# Author: ChatGPT
# Date: 2025
###########################################################

# ----------- PARAMETERS -----------------
INPUT_DIR="/home/zl436/palmer_scratch/bam/ATAC_final"
OUTPUT_DIR="/home/zl436/palmer_scratch/ATAC_analysis/fragment_lengths"

module load SAMtools

# Create output directory
mkdir -p $OUTPUT_DIR

# ----------- PROCESS EACH BAM -----------
echo "Extracting fragment length distributions..."

for BAM in $INPUT_DIR/*.bam; do
    SAMPLE=$(basename $BAM .bam)
    echo "Processing $SAMPLE..."

    # Extract TLEN field (column 9), keep only properly paired reads
    # Remove header, filter for non-zero, positive insert sizes
    samtools view -f 0x2 $BAM | awk '{if($9>0) print $9}' > $OUTPUT_DIR/${SAMPLE}_frag_lengths.txt

    echo "Saved fragment lengths to $OUTPUT_DIR/${SAMPLE}_frag_lengths.txt"
done

echo "✅ All samples processed!"
