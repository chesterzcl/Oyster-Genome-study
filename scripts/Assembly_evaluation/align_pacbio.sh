#!/bin/bash
#SBATCH --job-name=align_pb_reads
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G


module load minimap2
module load SAMtools

GENOME="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
READS="/home/zl436/ycga_work/oyster_genome/draft_genome/PacBio/hifi_reads/11S_pb.fastq.gz"
OUT_PREFIX="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/hifi_alignment"

if [ ! -f "${GENOME}.mmi" ]; then
    echo "Indexing genome..."
    minimap2 -d ${GENOME}.mmi ${GENOME}
fi

minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${GENOME}.mmi ${READS} | \
samtools view -@ 4 -bS - | samtools sort -@ 4 -o ${OUT_PREFIX}.bam


samtools index ${OUT_PREFIX}.bam

