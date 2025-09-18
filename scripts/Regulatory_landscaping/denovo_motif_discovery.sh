#!/bin/bash
#SBATCH --job-name=homer_motif
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=ALL

WORKDIR=/home/zl436/palmer_scratch/script/homer_dir
BED_FILE=consensus_gene.bed
GENOME_FASTA=primary_dedup_chr_masked_hp_sealed.fa
OUTPUT_DIR=gene_vs_concensus
GENOME_NAME=customOyster
SIF_IMAGE=homer.sif

cd ${WORKDIR}

# === Set up HOMER genome dir (host side) ===

GENOME_DIR=${WORKDIR}/.homer/${GENOME_NAME}
mkdir -p ${GENOME_DIR}
cp -u ${WORKDIR}/${GENOME_FASTA} ${GENOME_DIR}/${GENOME_NAME}.fa

# === Pull image if not exists ===
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "Pulling HOMER container..."
    apptainer build ${SIF_IMAGE} docker:quay.io/biocontainers/homer:5.1--pl5262h9948957_0
fi

# === Run HOMER inside Apptainer ===
apptainer exec \
  ${SIF_IMAGE} \
  bash -c "
    findMotifsGenome.pl \
      ${WORKDIR}/${BED_FILE} \
      ${WORKDIR}/${GENOME_FASTA} \
      ${WORKDIR}/${OUTPUT_DIR} \
      -bg consensus_annotated.bed \
      -len 8,10,12 \
      -size given \
      -p 4 
  "

