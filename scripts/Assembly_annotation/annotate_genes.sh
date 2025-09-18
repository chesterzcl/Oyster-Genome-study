#!/bin/bash
#SBATCH --job-name=annotate_gene
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G


# -------------------------------
# User-defined paths
DIR=/home/zl436/palmer_scratch/ref
PROTEIN_FASTA="${DIR}/gene_pred_mg/proteins.fa"
OUTDIR="${DIR}/gene_pred_mg/eggnog"
EMAPPER_DB_DIR="/home/zl436/palmer_scratch/script/emapper_dir/eggnog-mapper_mg/data"
CPU=16
CONTAINER="eggnog-mapper.sif"
# -------------------------------

apptainer pull $CONTAINER docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_2

# 2. Create output directory
mkdir -p "$OUTDIR"

mkdir -p "$EMAPPER_DB_DIR"


# 3. Download the eggNOG DB if not present
# Run the download script inside the container
apptainer exec \
  --bind "$EMAPPER_DB_DIR":/db \
  "$CONTAINER" \
  python /opt/eggnog-mapper/download_eggnog_database.py --data_dir /db

apptainer exec \
  --bind $(dirname "$PROTEIN_FASTA"):/input \
  --bind "$OUTDIR":/output \
  --bind "$EMAPPER_DB_DIR":/db \
  "$CONTAINER" \
  emapper.py \
    -i /input/$(basename "$PROTEIN_FASTA") \
    -o eggnog_annotation \
    -m diamond \
    --cpu $CPU \
    --output_dir /output \
    --data_dir /db \
    --override
