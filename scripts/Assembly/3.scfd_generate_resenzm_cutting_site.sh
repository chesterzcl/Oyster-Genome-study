#!/bin/bash
#SBATCH --job-name=HiC
#SBATCH --time=40:00:00
#SBATCH --partition=week
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=all

module load Java
module load GCC


apptainer pull juicer.sif docker://aidenlab/juicer:v2.0.1

dir=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp
fa_file=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa
tool_dir=/home/zl436/palmer_scratch/script/HiC_pipe_juicer

apptainer exec juicer.sif bwa index ${fa_file}

#apptainer exec juicer.sif find / -name generate_site_positions.py 2>/dev/null

apptainer exec juicer.sif python /opt/juicer/misc/generate_site_positions.py Arima final ${fa_file}

echo 'Restriction enzyme site file created'

awk 'BEGIN{OFS="\t"}{print $1, $NF}' final_Arima.txt > Arima_final.sizes
apptainer exec juicer.sif juicer.sh -t 12 --assembly -z ${fa_file} -p ${tool_dir}/final_Arima.sizes -y ${tool_dir}/final_Arima.txt

