#!/bin/bash
#SBATCH --job-name=rnaseq_pipeline
#SBATCH --time=23:59:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5G

module load Java
# Run the pipeline
./nextflow run main.nf -profile slurm -resume
