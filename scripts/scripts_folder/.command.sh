#!/bin/bash -ue
ls 11S_1_trimmed_grmv.fastq.gz 11S_2_trimmed_grmv.fastq.gz>sgs.fofn
	ls 11S_pb.fastq.gz primary.fa>hifi.fofn

cat << 'EOF' > run.cfg
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = ./primary.fa
genome_size = auto
workdir = ./01_rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa

[hifi_option]
hifi_fofn = ./hifi.fofn
hifi_options = -min_read_len 1k -max_depth 100
hifi_minimap2_options = -x map-pb
EOF


	nextPolish run.cfg
