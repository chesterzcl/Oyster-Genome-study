#!/bin/bash
#SBATCH --job-name=purge_dups
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL

module load minimap2
module load SAMtools
module load GCC

genome="/home/zl436/ycga_work/oyster_genome/draft_genome/Contigs/primary.fa"
hifi_reads="/home/zl436/ycga_work/oyster_genome/draft_genome/PacBio/hifi_reads/11S_pb.fastq.gz"  # HiFi reads
tool_dir=/vast/palmer/scratch/hoh/zl436/script/purge_dups/src
op_prefix=/vast/palmer/scratch/hoh/zl436/script/primary_dedup

minimap2 -xasm20 ${genome} ${hifi_reads} | gzip -c - > ${op_prefix}.paf.gz

${tool_dir}/pbcstat ${op_prefix}.paf.gz
${tool_dir}/calcuts /vast/palmer/scratch/hoh/zl436/script/PB.stat > ${op_prefix}.cutoffs 2>${op_prefix}.calcults.log

${tool_dir}/split_fa ${genome} > ${op_prefix}.split.fa
minimap2 -xasm5 -DP ${op_prefix}.split.fa ${op_prefix}.split.fa | gzip -c - > ${op_prefix}.split.self.paf.gz

${tool_dir}/purge_dups -2 -T ${op_prefix}.cutoffs -c /vast/palmer/scratch/hoh/zl436/script/PB.base.cov ${op_prefix}.split.self.paf.gz > ${op_prefix}.dups.bed 2> ${op_prefix}.purge_dups.log
${tool_dir}/get_seqs -e ${op_prefix}.dups.bed ${genome} 

