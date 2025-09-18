#!/bin/bash
#SBATCH --job-name=TOBIAS_TFBScan_Score
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

echo "==========================================="
echo "       TOBIAS pipeline: TFBScan + ScoreBed "
echo "       STARTED AT $(date)"
echo "==========================================="

# -------------------------------
# Load your environment
# -------------------------------
module load SAMtools
module load miniconda
conda activate seq_env

module load BEDTools

# -------------------------------
# USER CONFIGURATION
# -------------------------------
WORKDIR="/home/zl436/palmer_scratch/script/tobias_dir"
GENOME_FASTA="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
#CONSENSUS_BED="/home/zl436/palmer_scratch/script/tobias_dir/SE_peaks.bed"
CONSENSUS_BED="/home/zl436/palmer_scratch/bed/ATAC_final/consensus_peaks.bed"
FOOTPRINT_BW="${WORKDIR}/Footprints_all_peaks.bw"
MOTIFS_FILE="${WORKDIR}/JASPAR_core.txt"

# Create output directories
SCAN_OUTDIR="${WORKDIR}/TFBScan_out"
mkdir -p "$SCAN_OUTDIR"

SCORE_OUTDIR="${WORKDIR}/ScoreBed_out"
mkdir -p "$SCORE_OUTDIR"

# -------------------------------
# Step 1: Extract FASTA from peaks
# -------------------------------
echo "[Step 1] Extracting FASTA from BED peaks..."
PEAKS_FASTA="${SCAN_OUTDIR}/SE_peaks.fa"

bedtools getfasta \
    -fi "$GENOME_FASTA" \
    -bed "$CONSENSUS_BED" \
    -fo "$PEAKS_FASTA"

echo "FASTA extraction complete: $PEAKS_FASTA"


# -------------------------------
# Step 2: TFBScan to find motif sites
# -------------------------------
echo "[Step 2] Running TOBIAS TFBScan..."

TFBSCAN_BED="${SCAN_OUTDIR}/MotifMatches.bed"

TOBIAS TFBScan \
    --motifs "$MOTIFS_FILE" \
    --fasta "${GENOME_FASTA}" \
    --outfile "${TFBSCAN_BED}" \
    --regions "${CONSENSUS_BED}" \
    --cores 8

echo "TFBScan complete: $TFBSCAN_BED"


# -------------------------------
# Step 3: ScoreBed to assign footprint scores
# -------------------------------
echo "[Step 3] Running TOBIAS ScoreBed..."

SCORED_BED="${SCORE_OUTDIR}/ScoredMotifs.bed"

TOBIAS ScoreBed \
    --bed "$TFBSCAN_BED" \
    --bigwigs "$FOOTPRINT_BW" \
    --output "$SCORED_BED"

echo "ScoreBed complete: $SCORED_BED"

echo "==========================================="
echo "       TOBIAS pipeline COMPLETE!"
echo "       ENDED AT $(date)"
echo "==========================================="


