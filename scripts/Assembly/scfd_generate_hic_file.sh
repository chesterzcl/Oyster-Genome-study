#!/bin/bash
#SBATCH --job-name=HiC
#SBATCH --time=12:00:00
#SBATCH --partition=day
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=all

module load Java
module load GCC

apptainer pull juicertools.sif docker://adthrasher/juicertools:1.6.2

dir=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp
fa_file=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa
tool_dir=/home/zl436/palmer_scratch/script/HiC_pipe_juicer

cd ${tool_dir}
java -jar juicer_tools.2.20.00.jar pre --threads 4 \
-q 30 ${tool_dir}/aligned/merged_nodups.txt \
${tool_dir}/primary_dedup_chr_masked_hp_sealed.hic \
${tool_dir}/Arima_final.sizes \
BP 5000,10000,25000,50000,100000


