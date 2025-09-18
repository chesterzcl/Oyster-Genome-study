#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL


module load SeqKit

dir=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc
INPUT_FASTA=${dir}/primary_dedup.FINAL.fasta
OUTPUT_FASTA=${dir}/primary_dedup_chr.fa

# Run seqkit to filter
echo "Filtering contigs > mb from $INPUT_FASTA ..."
seqkit seq -m 20000000 "$INPUT_FASTA" > "$OUTPUT_FASTA"

echo "Done! Filtered genome saved to $OUTPUT_FASTA"

