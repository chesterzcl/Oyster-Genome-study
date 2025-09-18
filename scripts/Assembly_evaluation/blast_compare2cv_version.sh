#!/bin/bash
#SBATCH --job-name=blast_search
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2

# Load BLAST+ module (adjust based on your cluster setup)
module load BLAST+

# Define input and output
ASSEMBLY="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/primary_dedup_chr.fa"
QUERY="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/CV_sig.fa"
DB_NAME="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/primary_dedup_chr_blast_db"
OUT_FILE="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/blast_op.tsv"


makeblastdb -in ${ASSEMBLY} -dbtype nucl -out ${DB_NAME}

blastn -query ${QUERY} \
       -db ${DB_NAME} \
       -out ${OUT_FILE} \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
       -num_threads ${SLURM_CPUS_PER_TASK}
