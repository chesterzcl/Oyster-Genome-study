#!/bin/bash


TRF_DAT="$1"
FAI="$2"  # .fai index of your genome (chrom sizes)
PREFIX="satellite_analysis"

# Output files
TRF_BED="${PREFIX}_all_trf.bed"
SAT_BED="${PREFIX}_satellite_filtered.bed"
MOTIF_COUNTS="${PREFIX}_motif_counts.tsv"
SAT_STATS="${PREFIX}_satellite_summary.txt"

echo "Converting TRF output to BED format..."
awk -v OFS="\t" '
BEGIN { chrom = "" }
/^Sequence:/ { chrom = $2; next }
/^[0-9]/ && NF >= 14 {
    start = $1 - 1;
    end = $2;
    period = $3;
    copies = $4;
    stl_match = $6;
    indel = $7;
    score = $8;
    consensus = $14;
    stl_length = end - start;

    print chrom, start, end, consensus, score, ".", period, copies, stl_match, indel, stl_length;
}' "$TRF_DAT" > "$TRF_BED"

echo "Filtering satellite-like repeats (period>=5, copies>=10, length>=100bp, match>=85%)..."
awk -v OFS="\t" '$7 >= 5 && $8 >= 10 && $9 >= 85 && $11 >= 100' "$TRF_BED" > "$SAT_BED"

echo "Counting satellite bases and most common motifs..."
TOTAL_BP=$(awk '{sum += $3 - $2} END {print sum}' "$SAT_BED")
GENOME_SIZE=$(awk '{sum += $2} END {print sum}' "$FAI")
PCT=$(awk -v b="$TOTAL_BP" -v g="$GENOME_SIZE" 'BEGIN {printf "%.2f", (b/g)*100}')

echo -e "Total satellite bp:\t$TOTAL_BP" > "$SAT_STATS"
echo -e "Genome size (from .fai):\t$GENOME_SIZE" >> "$SAT_STATS"
echo -e "Percent satellite:\t$PCT%" >> "$SAT_STATS"

echo "Top 20 repeat motifs:"
cut -f4 "$SAT_BED" | sort | uniq -c | sort -nr | head -20 > "$MOTIF_COUNTS"
cat "$MOTIF_COUNTS"

echo "Satellite BED file ready: $SAT_BED"
echo "Summary file: $SAT_STATS"
echo "Motif frequencies: $MOTIF_COUNTS"

