#!/bin/bash
#SBATCH --job-name=gffread_extract
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00


apptainer pull gffread.sif docker://dceoy/gffread

# Set file paths
DIR=/home/zl436/palmer_scratch/ref
GENOME="${DIR}/Mgiga.fa"
GFF="${DIR}/Mgiga.gff"
OUTDIR="${DIR}/gene_pred_mg"

mkdir -p $OUTDIR

apptainer exec \
  --bind $GENOME:/data/genome.fa \
  --bind $GFF:/data/annotation.gff3 \
  --bind $OUTDIR:/output \
  gffread.sif \
  gffread /data/annotation.gff3 \
    -g /data/genome.fa \
    -y /output/proteins.fa \
    -x /output/cds.fa \
    -w /output/transcripts.fa

