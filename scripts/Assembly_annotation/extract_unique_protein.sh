#!/usr/bin/env bash

set -euo pipefail

############################################################
# USER CONFIGURATION
############################################################
INPUT="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/gene_pred/proteins.fa"
OUTPUT="/home/zl436/palmer_scratch/script/blast_dir/CV_final_protein_unique.fa"

############################################################
echo "Input : $INPUT"
echo "Output: $OUTPUT"
echo "Extracting longest isoforms per gene and sorting..."

############################################################
# Extract, pick longest isoform per gene, sort by gene ID
############################################################

awk '
    BEGIN { prefix = "" }
    /^>/ {
        if (seq_id != "") {
            gene = seq_id
            sub(/\..*/, "", gene)
            if (length(seq) > maxlen[gene]) {
                maxlen[gene] = length(seq)
                header[gene] = prefix gene
                sequence[gene] = seq
            }
        }
        seq_id = substr($0, 2)
        seq = ""
        next
    }
    {
        seq = seq $0
    }
    END {
        if (seq_id != "") {
            gene = seq_id
            sub(/\..*/, "", gene)
            if (length(seq) > maxlen[gene]) {
                maxlen[gene] = length(seq)
                header[gene] = prefix gene
                sequence[gene] = seq
            }
        }
        for (g in header) {
            print g "\t" header[g] "\t" sequence[g]
        }
    }
' "$INPUT" | sort -t$'\t' -k1,1V | awk -F'\t' '{print ">"$2"\n"$3}' > "$OUTPUT"

############################################################


