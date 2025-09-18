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
ATAC_PEAKS_BED="/home/zl436/palmer_scratch/bed/ATAC_final/consensus_peaks.bed"
OUTDIR="/home/zl436/palmer_scratch/script/tobias_dir"

mkdir -p ${OUTDIR}/TF_motifs

TOBIAS BINDetect \
    --motifs ${OUTDIR}/JASPAR_core.txt \
    --signals ${OUTDIR}/Footprints_all_peaks.bw \
    --peaks ${ATAC_PEAKS_BED} \
    --genome ${GENOME_FASTA} \
    --outdir ${OUTDIR}/TF_motifs \
    --cores $SLURM_CPUS_PER_TASK

