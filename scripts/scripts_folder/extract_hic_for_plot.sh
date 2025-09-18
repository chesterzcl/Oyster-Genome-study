#!/bin/bash
#SBATCH --job-name=hic_extract
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=ALL

# Define variables
SIF=/home/zl436/palmer_scratch/script/HiC_pipe_juicer/juicertools.sif
JAR=/mnt/jars/juicer_tools.2.20.00.jar
JAR_HOST=/home/zl436/palmer_scratch/script/HiC_pipe_juicer/juicer_tools.2.20.00.jar
HIC_FILE=/home/zl436/palmer_scratch/script/HiC_pipe_juicer/primary_dedup_chr_masked_hp_sealed.hic
RES=10000
OUT=ALL_${RES}_matrix.txt

CHROMS=(HiC_scaffold_1 HiC_scaffold_2 HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 \
        HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_8 HiC_scaffold_9 HiC_scaffold_10)

# === Initialize output file ===
> $OUT

# === Dump and combine by chromosome pairs ===
for chr1 in "${CHROMS[@]}"; do
  for chr2 in "${CHROMS[@]}"; do
    echo "Dumping: $chr1 vs $chr2"
    apptainer exec --bind $(dirname $JAR_HOST):/mnt/jars $SIF \
      java -Xmx4g -jar $JAR dump observed NONE $HIC_FILE $chr1 $chr2 BP $RES tmp.txt
    if [ -s tmp.txt ]; then
      awk -v c1=$chr1 -v c2=$chr2 '{print c1"\t"c2"\t"$1"\t"$2"\t"$3}' tmp.txt >> $OUT
    else
      echo "Warning: empty matrix for $chr1 vs $chr2"
    fi
    rm -f tmp.txt
  done
done


