#!/bin/bash
#SBATCH --job-name=SCFD
#SBATCH --time=48:00:00
#SBATCH --partition=week
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=all

module load Java
module load GCC

apptainer pull 3d-dna.sif docker://aakashsur/3d-dna

dir=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup
fa_file=${dir}/primary_dedup.fa
mnd_file=${dir}/merged_nodups.txt
tool_dir=/home/zl436/palmer_scratch/simu/tool/HiC_pipe_3ddna

apptainer exec 3d-dna.sif 3d-dna -r 0 ${fa_file} ${mnd_file}

