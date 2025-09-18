#!/bin/bash

# === USAGE EXAMPLE ===
# bash get_gc_content.sh genome.fa 100000 output_gc_content.tsv

# === INPUTS ===
FASTA="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
BIN_SIZE=200000
OUT_TSV="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/gc.tsv"

module load BEDTools
module load SAMtools

if [ -z "$FASTA" ] || [ -z "$BIN_SIZE" ] || [ -z "$OUT_TSV" ]; then
  echo "Usage: $0 <genome_fasta> <bin_size> <output_tsv>"
  exit 1
fi

# === 1. Make sure FASTA is indexed ===
if [ ! -f "${FASTA}.fai" ]; then
  echo "[INFO] Indexing FASTA..."
  samtools faidx "$FASTA"
fi

# === 2. Make BED windows ===
echo "[INFO] Creating bins of size ${BIN_SIZE}..."
bedtools makewindows -g "${FASTA}.fai" -w $BIN_SIZE > bins_${BIN_SIZE}.bed

# === 3. Compute GC content with bedtools nuc ===
echo "[INFO] Computing GC content..."
bedtools nuc -fi "$FASTA" -bed bins_${BIN_SIZE}.bed > bins_${BIN_SIZE}_nuc.tsv

# === 4. Extract relevant columns ===
# bedtools nuc output has:
# 1=chrom  2=start  3=end  ...  5=%GC
echo "[INFO] Extracting GC content..."
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3,$5}' bins_${BIN_SIZE}_nuc.tsv > "$OUT_TSV"

echo "[DONE] GC-content table written to $OUT_TSV"
