#!/bin/bash
#SBATCH --job-name=scan_motif
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# ===== USER PARAMETERS =====
WORKDIR=/home/zl436/palmer_scratch/script/homer_dir
GENOME_FASTA=primary_dedup_chr_masked.fa
MOTIF_FILE=homer_output_transposon/dna_transposon_top10.motifs            # HOMER-formatted motif file
GENOME_NAME=customOyster
OUTPUT_FILE=transposon_motif_hits.bed
SIF_IMAGE=homer.sif

cd ${WORKDIR}

# ===== Set up HOMER genome directory structure =====
GENOME_DIR=${WORKDIR}/.homer/${GENOME_NAME}
mkdir -p ${GENOME_DIR}
cp -u ${WORKDIR}/${GENOME_FASTA} ${GENOME_DIR}/${GENOME_NAME}.fa

# ===== Pull container if missing =====
if [ ! -f "${SIF_IMAGE}" ]; then
    echo "Pulling HOMER container..."
    apptainer build ${SIF_IMAGE} docker://quay.io/biocontainers/homer:5.1--pl5262h9948957_0
fi

# ===== Run motif scanning =====
apptainer exec \
  ${SIF_IMAGE} \
  bash -c "
    scanMotifGenomeWide.pl \
      ${WORKDIR}/${MOTIF_FILE} \
      ${WORKDIR}/${GENOME_FASTA} \
      -bed \
      -keepAll \
      > ${WORKDIR}/${OUTPUT_FILE}
  "

echo "Motif scanning done. Results written to: ${OUTPUT_FILE}"
