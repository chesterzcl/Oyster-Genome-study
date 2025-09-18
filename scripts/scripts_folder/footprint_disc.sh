#!/bin/bash
#SBATCH --job-name=TOBIAS_merged
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G


#######################################
module load SAMtools
module load miniconda
conda activate seq_env

#######################################
# USER CONFIGURATION
#######################################

INPUT_BAM_DIR="/home/zl436/palmer_scratch/bam/ATAC_final"
MERGED_BAM="merged.bam"
GENOME_FASTA="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
#ATAC_PEAKS_BED="/home/zl436/palmer_scratch/bed/ATAC_final/consensus_peaks.bed"
ATAC_PEAKS_BED="/home/zl436/palmer_scratch/script/tobias_dir/SE_peaks.bed"
OUTDIR="/home/zl436/palmer_scratch/script/tobias_dir"

# Make output dir
mkdir -p "$OUTDIR"

#######################################
# Step 1: Merge BAMs
#######################################
echo "Merging BAM files..."
samtools merge -@ $SLURM_CPUS_PER_TASK "$OUTDIR/$MERGED_BAM" "$INPUT_BAM_DIR"/*.bam
samtools index "$OUTDIR/$MERGED_BAM"

#######################################
# Step 2: ATACorrect
#######################################
echo "Running ATACorrect..."
TOBIAS ATACorrect \
    --bam "$OUTDIR/$MERGED_BAM" \
    --genome "$GENOME_FASTA" \
    --outdir "$OUTDIR" \
    --cores $SLURM_CPUS_PER_TASK \
    --peaks ${ATAC_PEAKS_BED}

# Resulting corrected signal:
CORRECTED_BW="$OUTDIR/${MERGED_BAM%.bam}_corrected.bw"

#######################################
# Step 3: FootprintScores
#######################################
echo "Running FootprintScores..."
TOBIAS FootprintScores \
    --signal "$CORRECTED_BW" \
    --regions "$ATAC_PEAKS_BED" \
    --output "$OUTDIR/Footprints_all_peaks.bw" \
    --cores $SLURM_CPUS_PER_TASK


