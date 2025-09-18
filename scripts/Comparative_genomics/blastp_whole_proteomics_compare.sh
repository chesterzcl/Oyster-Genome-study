#!/bin/bash
#SBATCH --job-name=blastp_mcscanx
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G

#####################################
# USER CONFIGURATION
#####################################

# Path to input FASTA files
QUERY_FASTA="/home/zl436/palmer_scratch/script/blast_dir/CV_final_vs_CV/CV_protein_unique.fa"
DB_FASTA="/home/zl436/palmer_scratch/script/blast_dir/CV_final_vs_CV/CV_final_protein_unique.fa"

# Output paths
DB_NAME="CV_final_db"
BLAST_OUT="CV_vs_CV_final.blast"

#####################################
# Load BLAST module (adapt to your cluster)
#####################################
module load BLAST+

#####################################
# Make BLAST database
#####################################
echo "Making BLASTP database..."
makeblastdb -in "$DB_FASTA" -dbtype prot -out "$DB_NAME"

#####################################
# Run BLASTP
#####################################
echo "Running BLASTP..."
blastp -query "$QUERY_FASTA" \
       -db "$DB_NAME" \
       -out "$BLAST_OUT" \
       -evalue 1e-5 \
       -outfmt 6 \
       -num_threads $SLURM_CPUS_PER_TASK


