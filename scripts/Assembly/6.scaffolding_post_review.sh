#!/bin/bash
#SBATCH --job-name=HiC
#SBATCH --time=48:00:00
#SBATCH --partition=week
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=all

module load Java
module load GCC

apptainer pull 3d-dna.sif docker://aakashsur/3d-dna

dir=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna
fa_file=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup/primary_dedup.fa
mnd_file=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup/merged_nodups.txt
tool_dir=/home/zl436/palmer_scratch/simu/tool/HiC_pipe_3ddna
asm_file=${dir}/primary_dedup_scfd/primary_dedup.rawchrom.review.assembly

apptainer exec 3d-dna.sif /root/3d-dna/run-asm-pipeline-post-review.sh --sort-output -r ${asm_file} ${fa_file} ${mnd_file}

