#!/bin/bash
#SBATCH --job-name=hifi_depth
#SBATCH --partition=day
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G

module load SAMtools

DIR="/home/zl436/palmer_scratch/bam/2nd_draft_final"
BAM="${DIR}/11S_st_dr_filtered.bam"
OUT="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/wgs_alignment_depth.txt"
BIN="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/wgs_depth_binned.txt"

echo "Running samtools depth on $BAM"

samtools depth -d 10000 "$BAM" > "$OUT"

echo "Depth calculation complete."

module load Python

source ~/my_python_module/bin/activate

python bin_depth.py \
  -i ${OUT} \
  -o ${BIN} \
  -b 100000
