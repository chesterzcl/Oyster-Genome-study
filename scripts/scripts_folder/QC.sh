#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL

module load FastQC

dir=/home/zl436/palmer_scratch/oyster_genome/PacBio/fail_reads

fastqc -t 8 ${dir}/oyster_pb.fastq.gz -o ${dir}/qc_results/
