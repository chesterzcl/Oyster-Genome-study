#!/bin/bash
#SBATCH --job-name=braker
#SBATCH --time=23:00:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=ALL

apptainer pull braker3.sif docker://teambraker/braker3:latest 

export AUGUSTUS_CONFIG_PATH=$PWD/augustus_config
mkdir -p $AUGUSTUS_CONFIG_PATH

input_dir=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc
genome=${input_dir}/primary_dedup_chr_masked.fa
bam_file=${input_dir}/rna_merged.bam
op_dir=${input_dir}/gene_pred_2

mkdir -p ${op_dir}

apptainer exec \
  --bind $PWD:/data \
  --bind $AUGUSTUS_CONFIG_PATH:/config \
  --bind ${input_dir}:${input_dir} \
  braker3.sif \
  bash -c "
    export AUGUSTUS_CONFIG_PATH=/config
    export GENEMARK_PATH=/opt/ETP/bin
    export AUGUSTUS_BIN_PATH=/opt/Augustus/bin
    export AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts
    export PYTHON3_PATH=/opt/conda/bin

    braker.pl \
      --genome=${genome} \
      --bam=${bam_file} \
      --softmasking \
      --threads=16 \
      --species=oyster \
      --workingdir=${op_dir} \
      --gff3 \
      --useexisting
  "


