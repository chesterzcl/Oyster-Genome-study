#!/bin/bash
#SBATCH --job-name=hic_stage1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G

module load miniconda
conda activate seq_env

# Resolution parameter
res=1000

# Base directory
BASEDIR="/home/zl436/palmer_scratch/script/HiC_explorer_dir"

# Input and output
INPUT_COOL="${BASEDIR}/ALL_${res}.cool"
OUTDIR="${BASEDIR}/analysis_${res}"
mkdir -p "$OUTDIR"

# Convert cool to h5
OUTPUT_H5="${OUTDIR}/ALL_${res}.h5"
echo "Converting ${INPUT_COOL} -> ${OUTPUT_H5}"
hicConvertFormat --matrices "$INPUT_COOL" --outFileName "$OUTPUT_H5" --inputFormat cool --outputFormat h5

# Make diagnostic plot
echo "Generating QC diagnostic plot..."
hicCorrectMatrix diagnostic_plot -m "$OUTPUT_H5" -o "${OUTDIR}/ALL_${res}_QC.png"

echo "Stage 1 complete!"
echo "Please review the QC plot at: ${OUTDIR}/ALL_${res}_QC.png"

