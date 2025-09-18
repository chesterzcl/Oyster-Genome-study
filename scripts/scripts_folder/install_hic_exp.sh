#!/bin/bash
#SBATCH --job-name=hic2cool
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

module load miniconda

conda activate seq_env

conda install hicexplorer -c bioconda -c conda-forge
#conda install -c bioconda hic2cool

