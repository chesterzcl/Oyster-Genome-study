#!/bin/bash
#SBATCH --job-name=satellite_detect
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

# ---------- User-configurable ----------
GENOME="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
TRF_PARAMS="2 7 7 80 10 50 2000"
OUTDIR="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/stlt_analysis"
MIN_PERIOD=10
MIN_LENGTH=1000

module load TRF


# Setup
mkdir -p "$OUTDIR"
cd "$OUTDIR"

BASENAME=${GENOME}
TRF_DAT="${BASENAME}.${TRF_PARAMS// /\.}.dat"
TRF_BED="${BASENAME}_trf_all.bed"
SAT_BED="${BASENAME}_satellites.bed"

# Step 1: Run TRF
echo "Running Tandem Repeat Finder..."
trf "$GENOME" $TRF_PARAMS -d -h

# Step 2: Convert TRF .dat output to BED with period and length
echo "Parsing TRF output to BED format..."
awk -v OFS="\t" '!/^#/ && NF >= 14 {
    chrom=$14;
    start=$1-1;
    end=$2;
    score=$3;
    period=$4;
    trf_length=end-start;
    strand=($13==1 ? "+" : "-");
    print chrom, start, end, "TRF_"NR, score, strand, period,trf_length
}' "$TRF_DAT" > "$TRF_BED"

# Step 3: Filter candidate satellites
echo "Filtering candidate satellite regions..."
awk -v OFS="\t" -v p="$MIN_PERIOD" -v l="$MIN_ARRAY_LEN" '$7 >= p && $8 >= l' "$TRF_BED" > "$SAT_BED"

