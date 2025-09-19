#!/bin/bash

INFILE="primary_dedup_chr_masked_hp_sealed.fa.out"
OUTFILE="all_repeats.bed"

echo "Converting $INFILE to BED with class/family -> $OUTFILE"

awk 'BEGIN{OFS="\t"} NR<=2 {next} 
     {
       strand = ($9=="C" ? "-" : $9);
       print $5, $6-1, $7, $10, ".", strand, $11
     }' "$INFILE" > "$OUTFILE"

echo "Done. BED file written to $OUTFILE"

