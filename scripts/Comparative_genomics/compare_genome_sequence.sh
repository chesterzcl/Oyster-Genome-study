#!/bin/bash
#SBATCH --job-name=synteny_align
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=week
#SBATCH --mail-type=ALL

module load miniconda
source ~/.bashrc
conda activate seq_env

# Paths to genomes
REF=/home/zl436/palmer_scratch/ref/CV_clean.fa
REF=/home/zl436/palmer_scratch/ref/Mgiga_clean_fixed.fa
QUERY=/home/zl436/palmer_scratch/ref/CV_new_hr_clean.fa


# Output prefix
PREFIX=mg_vs_cv_new_all

# === Generate genome info file for plotsr ===
GENOME_INFO=${PREFIX}.txt
echo -e "#file\tname\ttags" > $GENOME_INFO
echo -e "$REF\tC.virginica_old\tlw:1.5" >> $GENOME_INFO
echo -e "$QUERY\tC.virginica_new\tlw:1.5" >> $GENOME_INFO


# Containers
MUMMER_IMG=mummer4_4.0.1.sif
SYRI_IMG=syri_1.7.0.sif

# Pull containers if not already present
if [ ! -f $MUMMER_IMG ]; then
    apptainer pull $MUMMER_IMG docker://quay.io/biocontainers/mummer4:4.0.1--pl5321h9948957_0
fi


# === Step 1: Whole-genome alignment ===
apptainer exec $MUMMER_IMG \
    nucmer --maxmatch -c 100 -b 500 -l 50 -t 16 $REF $QUERY -p $PREFIX


# === Step 2: Filter alignments ===
apptainer exec $MUMMER_IMG \
    delta-filter -m -i 90 -l 1000 ${PREFIX}.delta > ${PREFIX}.filtered.delta


# === Step 3: Convert delta to coords ===
apptainer exec $MUMMER_IMG \
    show-coords -THrd ${PREFIX}.filtered.delta > ${PREFIX}.filtered.coords


# === Step 4: Run SyRI using coords/delta ===
syri -c ${PREFIX}.filtered.coords \
-d ${PREFIX}.filtered.delta -f \
-r $REF -q $QUERY --prefix ${PREFIX}_ \
--no-chrmatch

echo "SyRI completed. Run plotsr separately using conda or native Python:"

plotsr --sr ${PREFIX}_syri.out \
--genomes $GENOME_INFO -H 8 -W 5 -o ${PREFIX}_plot.png --chrord chrorder.txt

conda deactivate

