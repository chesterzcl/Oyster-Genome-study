#!/usr/bin/env bash

INPUT_CV="CV_final.bed"
INPUT_MG="Mgiga.bed"
OUTPUT="combined.bed"

echo "✅ Extracting columns 1, 2, 4 and combining..."

(awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4, $2}' "$INPUT_CV"
 awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4, $2}' "$INPUT_MG") > "$OUTPUT"

echo "✅ Combined file written to $OUTPUT"


