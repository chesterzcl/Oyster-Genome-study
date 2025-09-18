params.data_dir="/home/zl436/palmer_scratch/bam/ATAC_final"
params.bams="${params.data_dir}/*_st_filtered.bam"
params.bam_idx="${params.data_dir}/*_st_filtered.bam.bai"
params.shifted_bams="${params.data_dir}/*_st_filtered_shifted.bam"


process Tn5_shift{
	container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'

	cpus 1

	input:
	tuple val(replicateID),path(bam_file),path(bam_index)

	output:
	tuple val(replicateID),path("${replicateID}_shifted.bam")

	script:
	"""
	alignmentSieve --bam ${bam_file} --outFile ${replicateID}_shifted.bam --ATACshift --numberOfProcessors 1
	"""

}



process ATAC_bam_QC{
	container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

	cpus 1

	input:
	tuple val(replicateID),path(bam_file),path(idx_file)

	output:
	tuple val(replicateID),path("${bam_file.baseName}_reads_sum.txt")

	script:
	"""
	samtools flagstat ${bam_file} > ${bam_file.baseName}_reads_sum.txt
	"""



}

process ATAC_bam_count_depth{

	container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
	
	cpus 1

        input:
        tuple val(replicateID),path(bam_file),path(idx_file)

        output:
        tuple val(replicateID),path("${bam_file.baseName}_depth.txt")

	script:
	"""
	samtools depth -a ${bam_file}>${bam_file.baseName}_depth.txt
	"""



}

process ATAC_peak_calling{
	container 'quay.io/biocontainers/macs2:2.1.1--r3.2.2_0'	

	cpus 1

	input:
	tuple val(replicateID),path(bam_file),path(idx_file)

	output:
	tuple val(replicateID),path("${bam_file.baseName}_peak_peaks.narrowPeak")

	script:
	"""
	macs2 callpeak -t ${bam_file} -f BAMPE --nomodel --shift -100 --extsize 200 -q 0.01 -n ${bam_file.baseName}_peak -g 53000000 --keep-dup all
	"""

}

process sort_index_bam{
	
        container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

        cpus 1

        input:
        tuple val(replicateID),path(bam_file)

        output:
        tuple val(replicateID),path("${replicateID}_shifted_mtrm_st.bam"),path("${replicateID}_shifted_mtrm_st.bam.bai")

        script:
        """
        samtools view -h ${bam_file} | grep -v 'NC_007175.2' | samtools view -bS -> ${bam_file.baseName}_mtrm.bam
	samtools sort -o ${bam_file.baseName}_mtrm_st.bam ${bam_file.baseName}_mtrm.bam
	samtools index ${bam_file.baseName}_mtrm_st.bam
        """

}

process count_per_peak_reads{
	container 'quay.io/biocontainers/bedtools:2.27.1--h077b44d_9'

	cpus 1

	input:
	tuple val(replicateID),path(bam),path(bam_idx),path(bed)
	
	output:
	tuple val(replicateID),path("${bam.baseName}_read_count.bed")
	
	script:
	"""
	bedtools multicov -bams ${bam} -bed ${bed} > ${bam.baseName}_read_count.bed
	"""

}


workflow{

	bams_ch=Channel.fromPath(params.bams,checkIfExists: true).map{file->tuple(file.baseName,file)}
	bams_ch.view()

	bam_idx_ch=Channel.fromPath(params.bam_idx).map{bai_file->tuple(bai_file.baseName.replace('.bam', ''),bai_file)}
	bam_idx_ch.view()

	bams_ch.join(bam_idx_ch).set{bam_with_idx_ch}
	bam_with_idx_ch.view()
		
	shifted_bam_ch=Tn5_shift(bam_with_idx_ch)
	shifted_bam_ch.view()

	shifted_bam_set_ch=shifted_bam_ch|sort_index_bam
	shifted_bam_set_ch.view()

        peak_files_ch=shifted_bam_set_ch|ATAC_peak_calling
        peak_files_ch.view()

	depth_ch=shifted_bam_set_ch|ATAC_bam_count_depth
	depth_ch.view()	

	sum_stats_ch=shifted_bam_set_ch|ATAC_bam_QC
	sum_stats_ch.view()	

	read_peak_ch=shifted_bam_set_ch.join(peak_files_ch)|count_per_peak_reads
	read_peak_ch.view()	

}
