#!/bin/bash
#SBATCH --job-name=repeatmodeler
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=week
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=16

# Paths
DIR=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp
GENOME="${DIR}/primary_dedup_chr_masked_hp_sealed.fa"
DBNAME="${DIR}/final_genome_db"
IMAGE="/home/zl436/palmer_scratch/script/repmdl_dir/tetools.sif"
THREADS=16


apptainer pull ${IMAGE} docker://dfam/tetools

# Create a working directory for RepeatModeler
mkdir -p repeatmodeler_output_final
cd repeatmodeler_output_final
cp $GENOME .

GENOME_LOCAL=$(basename "$GENOME")

# Step 1: Build RepeatModeler database
apptainer exec $IMAGE BuildDatabase -name $DBNAME $GENOME_LOCAL

# Step 2: Run RepeatModeler with LTRStruct
apptainer exec $IMAGE RepeatModeler -database $DBNAME -threads $THREADS -LTRStruct

# Step 3: Run RepeatMasker using the custom repeat library
# This assumes output library is ${DBNAME}-families.fa
apptainer exec $IMAGE RepeatMasker -pa $THREADS -lib ${DBNAME}-families.fa -xsmall $GENOME_LOCAL


