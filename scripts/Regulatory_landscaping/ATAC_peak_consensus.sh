#!/bin/bash
#SBATCH --job-name=consensus_peaks
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

cd /home/zl436/palmer_scratch/bed

# ======= PARAMETERS =======
MERGE_DIST=100       # Max distance between peaks to merge (in bp)
MIN_SUPPORT=3        # Minimum samples supporting a consensus peak
PEAK_SUFFIX="_peak_peaks.narrowPeak"

# ======= LOAD MODULES =======
module load BEDtools

# ======= STEP 1: Merge all peaks to get consensus regions =======
cat *${PEAK_SUFFIX} | sort -k1,1 -k2,2n > all_peaks.sorted.bed

bedtools merge -i all_peaks.sorted.bed -d $MERGE_DIST > consensus_peaks.bed

# ======= STEP 2: Prepare list of narrowPeak files =======
ls *${PEAK_SUFFIX} > peak_list.txt

# ======= STEP 3: Count sample support per consensus peak =======
bedtools multiinter -i $(cat peak_list.txt) > consensus_support.raw.bed

# Output format: chrom, start, end, N, list of files supporting

# ======= STEP 4: Filter by minimum support =======
awk -v min_supp=$MIN_SUPPORT 'BEGIN{OFS="\t"} {if ($4>=min_supp) print $1,$2,$3,$4}' consensus_support.raw.bed > consensus_peaks_min${MIN_SUPPORT}.bed

echo "Done! Output: consensus_peaks_min${MIN_SUPPORT}.bed"

# ======= STEP 5: Merge consecutive intervals =======
awk 'BEGIN{OFS="\t"} 
{
  if (prev && $1==chr && $2==end) {
    end=$3
  } else {
    if (prev) print chr, start, end, val
    chr=$1; start=$2; end=$3; val=$4; prev=1
  }
} 
END{print chr, start, end, val}' consensus_peaks_min${MIN_SUPPORT}.bed > consensus_peaks_min${MIN_SUPPORT}_merged.bed

echo "Done! Final merged output: consensus_peaks_min${MIN_SUPPORT}_merged.bed"





