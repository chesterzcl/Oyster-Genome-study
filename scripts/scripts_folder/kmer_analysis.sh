#!/bin/bash
#SBATCH --job-name=kmer
#SBATCH --time=3:00:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL

apptainer pull meryl.sif docker://quay.io/biocontainers/meryl:1.4.1--h9948957_2
apptainer pull merqury.sif docker://danylmb/merqury:1.4.1-build2

genome_dir=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc
genome=${genome_dir}/primary_dedup_chr.fa
hifi_reads=/home/zl436/palmer_scratch/oyster_genome/PacBio/hifi_reads/11S_pb.fastq.gz
kmer_dir=/home/zl436/palmer_scratch/kmer_db
kmer_db_wgs=${kmer_dir}/wgs_kmer.meryl
kmer_db_hifi=${kmer_dir}/hifi_kmer.meryl
kmer_db_asm=${kmer_dir}/draft2_kmer.meryl
kmer_db_draft2_hifi=${kmer_dir}/draft2_hifi_kmer.meryl
kmer_db_draft2_wgs=${kmer_dir}/draft2_wgs_kmer.meryl

#apptainer exec meryl.sif meryl statistics $kmer_db_asm>kmer_draft2_stats.txt

#apptainer exec meryl.sif meryl intersect $kmer_db_wgs $kmer_db_asm output $kmer_db_draft2_wgs
#apptainer exec meryl.sif meryl statistics $kmer_db_draft2_wgs>kmer_draft2_wgs_stats.txt

#apptainer exec meryl.sif meryl count k=21 output ${kmer_db_hifi} ${hifi_reads}
#apptainer exec meryl.sif meryl statistics $kmer_db_hifi>kmer_hifi_stats.txt

#apptainer exec meryl.sif meryl intersect $kmer_db_hifi $kmer_db_asm output $kmer_db_draft2_hifi
#apptainer exec meryl.sif meryl statistics $kmer_db_draft2_hifi>kmer_draft2_hifi_stats.txt

# Run Merqury for HiFi-based evaluation
apptainer exec merqury.sif merqury.sh $kmer_db_hifi $genome_dir/primary_dedup_chr.fa merqury_draft2_hifi
apptainer exec merqury.sif merqury.sh $kmer_db_wgs $genome_dir/primary_dedup_chr.fa merqury_draft2_wgs
