params.dir="/home/zl436/palmer_scratch/oyster_genome"
params.reads_hifi = "${params.dir}/PacBio/hifi_reads/11S_pb.fastq.gz"
params.reads_illumina = "${params.dir}/WGS/11S_*.fastq.gz"
params.reads_hic = "${params.dir}/HiC/11S_*_HiC.fastq.gz"


process pacbio_denovo{
	container 'quay.io/biocontainers/hifiasm:0.24.0--h5ca1c30_0'
	
	cpus 16
	
	input:
	path hifi_reads

	output:
	path "contigs*"

	script:
	"""
	hifiasm -o contigs -t 16 ${hifi_reads}
	"""

}

workflow{
	pb_reads=Channel.fromPath("${params.reads_hifi}")
	pacbio_denovo(pb_reads)
	pacbio_denovo.out.view()

}

