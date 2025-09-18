#!/bin/bash
#SBATCH --job-name=OrthoFinder
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

module load miniconda
# Load OrthoFinder module or activate your conda environment
source ~/.bashrc
conda activate orthofinder_env

# Step 1: Clean the protein FASTA files to remove '.' characters
CLEAN_DIR="cv_mg"

# Step 2: Run OrthoFinder
orthofinder -f $CLEAN_DIR -t 16 -a 16
