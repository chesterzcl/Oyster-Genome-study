#!/bin/bash

#SBATCH --time=12:00:00                 # Time limit (hh:mm:ss)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=32G                       # Memory per node

# Define input variables
BAM_DIR="/home/zl436/palmer_scratch/ATAC_L/bam"        # Directory containing BAM files
PEAK_DIR="/home/zl436/palmer_scratch/ATAC_L/peak"  # Directory for MACS2 output
BED_DIR="/home/zl436/palmer_scratch/ATAC_L/bed"
MIN_PEAK_SIZE=30           # Minimum peak size to retain
CONSENSUS_FILE="consensus_sorted.bed"

module load BEDTools

cat $PEAK_DIR/*_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > ${BED_DIR}/consensus_peaks.bed

awk -v min_size=$MIN_PEAK_SIZE '{if ($3-$2 > min_size) print $0}' ${BED_DIR}/consensus_peaks.bed > ${BED_DIR}/consensus_filtered.bed

sort -k1,1 -k2,2n ${BED_DIR}/consensus_filtered.bed > ${BED_DIR}/$CONSENSUS_FILE

bedtools multicov -bams $BAM_DIR/*.bam -bed ${BED_DIR}/$CONSENSUS_FILE > ${BED_DIR}/counts.txt

