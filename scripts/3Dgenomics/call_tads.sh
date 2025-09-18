#!/bin/bash
#SBATCH --job-name=hic_stage2
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G

module load miniconda
conda activate seq_env

# Resolution parameter
res=25000

# Choose your threshold (update after reviewing QC plot)
LOWER=-2.4
UPPER=4

# Base directory
BASEDIR="/home/zl436/palmer_scratch/script/HiC_explorer_dir"
OUTDIR="${BASEDIR}/analysis_${res}"

# Input/output
INPUT_H5="${OUTDIR}/ALL_${res}.h5"
CORRECTED_H5="${OUTDIR}/ALL_${res}_corrected.h5"

echo "Using thresholds: ${LOWER} to ${UPPER}"
echo "Correcting matrix..."
hicCorrectMatrix correct -m "$INPUT_H5" --filterThreshold $LOWER $UPPER -o "$CORRECTED_H5"

echo "Plotting corrected heatmap..."
hicPlotMatrix -m "$CORRECTED_H5" -o "${OUTDIR}/ALL_${res}_hm.png" --log1p

echo "Calling TADs..."
hicFindTADs -m "$CORRECTED_H5" \
  --outPrefix "${OUTDIR}/ALL_${res}_TADs" \
  --correctForMultipleTesting fdr \
  --numberOfProcessors 4 \
  --minDepth 100000 \
  --maxDepth 500000 \
  --step 50000 \
  --delta 0.005


echo "Stage 2 complete!"



