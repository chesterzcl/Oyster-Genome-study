#!/bin/bash
#SBATCH --job-name=bin_depth
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:50:00


module load Python

source ~/my_python_module/bin/activate

python bin_depth.py \
  -i /home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/hifi_depth.txt \
  -o /home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/hifi_bin_depth.tsv \
  -b 100000


