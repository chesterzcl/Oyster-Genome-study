#!/bin/bash

# === USAGE ===
# This script splits a genome FASTA into:
#   - Primary (haploid-like) assembly with retained regions
#   - Alternate haplotigs with removed regions

# --- INPUTS ---
FASTA=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/primary_dedup_chr_masked.fa
REMOVE_BED=/home/zl436/palmer_scratch/bed/haplotig_final.bed
OUTPUT_DIR=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/

PRIMARY_OUTPUT=${OUTPUT_DIR}/primary_dedup_chr_masked_hp.fa
ALTERNATE_OUTPUT=${OUTPUT_DIR}/alternate_haplotigs.fa

# --- Check index ---
if [ ! -f "${FASTA}.fai" ]; then
    echo "Index not found. Creating FASTA index..."
    samtools faidx "$FASTA"
fi

module load BEDTools

# --- Make sure output dir exists ---
mkdir -p "$OUTPUT_DIR"

# --- Step 1: Get complement intervals (regions to keep in primary) ---
echo "Computing retained regions..."
#bedtools complement -i "$REMOVE_BED" -g "${FASTA}.fai" > "${OUTPUT_DIR}/retained_regions.bed"

# --- Step 2: Extract primary (retained) regions ---
echo "Extracting primary (retained) sequences..."
#bedtools getfasta -fi "$FASTA" -bed "${OUTPUT_DIR}/retained_regions.bed" -fo "$PRIMARY_OUTPUT"

# --- Step 3: Extract alternate (removed) regions as haplotigs ---
echo "Extracting alternate (haplotig) sequences..."
bedtools getfasta -fi "$FASTA" -bed "$REMOVE_BED" -fo "$ALTERNATE_OUTPUT"

echo "Done."
echo "Primary assembly written to: $PRIMARY_OUTPUT"
echo "Alternate haplotigs written to: $ALTERNATE_OUTPUT"

