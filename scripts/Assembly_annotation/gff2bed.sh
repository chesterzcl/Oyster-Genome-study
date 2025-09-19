#!/usr/bin/env bash

INPUT_GFF="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/gene_pred/braker.gff3"
OUTPUT_BED="CV_final.bed"
WINDOW=2000

awk -F'\t' 'BEGIN {OFS="\t"}
  $3 == "gene" {
    chrom=$1
    sub(/^HiC_scaffold_/, "", chrom)
    chrnum=chrom+0
    start=$4
    end=$5
    split($9, a, "[=;]")
    id=a[2]
    sub(/^C\.Virginica_new\|/, "", id)
    print chrnum, start, "cf"chrnum, start, end, id
  }
' "$INPUT_GFF" | sort -n -k1,1 -k2,2 | cut -f3- > "$OUTPUT_BED"


