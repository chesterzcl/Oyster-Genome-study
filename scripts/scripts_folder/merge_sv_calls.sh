#!/bin/bash
#SBATCH --job-name=sv_merge
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Load apptainer module if needed (uncomment if your cluster uses modules)
# module load apptainer

# === Settings ===
SURVIVOR_IMAGE="survivor_1.0.7.sif"
CONTAINER_SOURCE="docker://quay.io/biocontainers/survivor:1.0.7--h077b44d_6"
SAMPLES_LIST="all_samples.txt"
DECOMP_DIR="decompressed_vcfs"
VCF_LIST="clean_vcf_list.txt"
MERGED_VCF="merged_sv.vcf"

# === Step 1: Pull SURVIVOR container if not present ===
if [[ ! -f "$SURVIVOR_IMAGE" ]]; then
    echo "Pulling SURVIVOR container..."
    apptainer pull "$SURVIVOR_IMAGE" "$CONTAINER_SOURCE"
fi

# === Step 2: Decompress each VCF ===
echo "Decompressing VCF files into: $DECOMP_DIR"
mkdir -p "$DECOMP_DIR"

while read vcf; do
    sample=$(basename "$(dirname "$(dirname "$(dirname "$vcf")")")")  # e.g. 10L_st_dr_filtered_manta_run
    out="$DECOMP_DIR/${sample}.vcf"
    zcat "$vcf" > "$out"
done < "$SAMPLES_LIST"

# === Step 3: Generate SURVIVOR input list ===
find "$DECOMP_DIR" -name "*.vcf" | sort > "$VCF_LIST"
echo "Generated list: $VCF_LIST"
cat "$VCF_LIST"

# === Step 4: Run SURVIVOR merge ===
echo "Running SURVIVOR merge..."
apptainer exec "$SURVIVOR_IMAGE" \
  SURVIVOR merge "$VCF_LIST" 1000 2 1 1 0 30 "$MERGED_VCF"

