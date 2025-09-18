#!/bin/bash
#SBATCH --job-name=HiC
#SBATCH --time=1:00:00
#SBATCH --partition=day
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=all


module load miniconda

conda activate seq_env

python hic2cool.py
