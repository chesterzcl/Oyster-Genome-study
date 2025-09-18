#!/bin/bash
#SBATCH --job-name=bam2fastq
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G


#apptainer pull bam2fastq.sif docker://quay.io/biocontainers/pbtk:3.5.0--h9ee0642_0

dir=/home/zl436/palmer_scratch/oyster_genome/PacBio/hifi_reads
bam_file="${dir}/11S_pb.bam"
out_prefix="${dir}/11S_pb_pbtk"

echo ${bam_file}
echo ${out_prefix}

/home/zl436/palmer_scratch/oyster_genome/pbtk/bam2fastq ${bam_file}

#apptainer exec bam2fastq.sif bam2fastq -o ${output_prefix} ${bam_file}
