params.genome_dir="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp"
params.genome="primary_dedup_chr_masked_hp_sealed.fa"
params.data_dir="/home/zl436/ycga_work/oyster_genome/ATAC/fastq"
params.reads="${params.data_dir}/*_ATAC_{1,2}.fastq.gz"


process prepare_genome_samtools{

        container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

        input:
        path genome

	output:
	path "${genome}.fai"

        script:
        """
        samtools faidx ${genome}
        """
}

process prepare_genome_bwa{

        container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

        input:
        path genome

	output:
	path "${genome.baseName}.*"

        script:
        """
        bwa index ${genome} -p ${genome.baseName}
        """
}

process prepare_genome_picard {

	container 'quay.io/biocontainers/picard:1.141--hdfd78af_6'
	
	input:
	path genome

	output:
	path "${genome.baseName}.dict"
	
	script:
	"""
	picard CreateSequenceDictionary R=${genome} O=${genome.baseName}.dict
	"""
}


process reads_quality_control{

	container 'quay.io/biocontainers/fastqc:0.11.9--0'
	
	cpus 5

	input:
	tuple val(replicateID),path(reads)
	
	output:
	tuple val(replicateID),path("${replicateID}_ATAC_1_fastqc.zip"),path("${replicateID}_ATAC_2_fastqc.zip")

	script:
	"""
	fastqc -t ${task.cpus} ${reads}
	"""
}

process reads_quality_control_post{

        container 'quay.io/biocontainers/fastqc:0.11.9--0'

        cpus 5

        input:
        tuple val(replicateID),path(reads)

        output:
        tuple val(replicateID),path("${replicateID}_1_trimmed_grmv_fastqc.zip"),path("${replicateID}_2_trimmed_grmv_fastqc.zip")

        script:
        """
        fastqc -t ${task.cpus} ${reads}
        """
}

process get_trimming_parameters{
	
	container 'chesterli214/gsexp_amd64_v0.6'

	input:
	tuple val(replicateID),path(QC_results_1),path(QC_results_2)

	output:
	tuple val(replicateID),env(trimming_params)

	script:
	"""
	unzip ${QC_results_1}
	unzip ${QC_results_2}
	qual1=\$(Get_trimming_parameter_qual ${QC_results_1.baseName}/fastqc_data.txt)
	qual2=\$(Get_trimming_parameter_qual ${QC_results_1.baseName}/fastqc_data.txt)
	trimming_params="ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:25 MINLEN:35"
	"""
}

process quality_trimming{

	container 'quay.io/biocontainers/trimmomatic:0.36--3'
	
	cpus 5

	input:
	tuple val(replicateID),path(reads),val(trimming_params)

	output:
	tuple val(replicateID),path("${replicateID}_*_trimmed.fastq.gz")	

	script:
	"""
	trimmomatic PE -threads ${task.cpus} ${reads} ${replicateID}_1_trimmed.fastq.gz ${replicateID}_1_unpaired.fastq.gz ${replicateID}_2_trimmed.fastq.gz ${replicateID}_2_unpaired.fastq.gz ${trimming_params}
	"""
}

process remove_polyG{
	container 'quay.io/biocontainers/fastp:0.24.0--heae3180_1'

	cpus 4

	input:
	tuple val(replicateID),path(reads)

	output:
        tuple val(replicateID),path("${replicateID}_*_trimmed_grmv.fastq.gz")


	script:
	"""
	fastp -i ${replicateID}_1_trimmed.fastq.gz -I ${replicateID}_2_trimmed.fastq.gz -o ${replicateID}_1_trimmed_grmv.fastq.gz -O ${replicateID}_2_trimmed_grmv.fastq.gz --trim_poly_g -l 35 --thread 4
	"""
}

process reads_mapping{
	container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

	cpus 5

	input:
	path genome_index
	tuple val(replicateID),path(trimmed_reads)

	output:
	tuple val(replicateID),path("${replicateID}.sam")

	script:
	"""
	bwa mem -t ${task.cpus} -R '@RG\\tID:${replicateID}\\tSM:${replicateID}.SM\\tLB:${replicateID}.LB' -o ${replicateID}.sam ${genome_index[0].baseName} ${trimmed_reads}
	rm ${trimmed_reads}
	"""
}

process post_alignment_processing{
	container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

	cpus 5

	input:
	tuple val(replicateID),path(samfile)

	output:
	tuple val(replicateID),path("${replicateID}_st_filtered.bam"),path("${replicateID}_st_filtered.bam.bai")

	script:
	"""
	samtools fixmate -@ ${task.cpus} -O BAM ${replicateID}.sam ${replicateID}.bam
	rm ${replicateID}.sam
	samtools sort -@ ${task.cpus} -O BAM -o ${replicateID}_sorted.bam ${replicateID}.bam
	rm ${replicateID}.bam
	samtools view -@ ${task.cpus} -O BAM -q 30 -o ${replicateID}_st_filtered.bam ${replicateID}_sorted.bam
	rm ${replicateID}_sorted.bam
	samtools index ${replicateID}_st_filtered.bam 
	"""
}



workflow{
        
	//***************************************************************
	//Generating genome indices

        prepare_genome_samtools("${params.genome_dir}/${params.genome}")

	prepare_genome_bwa("${params.genome_dir}/${params.genome}")

	//***************************************************************
	//Prepare input data

	reads_ch=Channel.fromFilePairs(params.reads,checkIfExists: true)
	reads_ch.view()

	//***************************************************************
	//Preprocessing

	reads_quality_control(reads_ch)
	
	get_trimming_parameters(reads_quality_control.out)
	get_trimming_parameters.out.view()

	reads_ch=Channel.fromFilePairs(params.reads,checkIfExists: true)
	
	reads_ch.join(get_trimming_parameters.out).set{QC_ready_reads}
	
	QC_ready_reads.view()
	quality_trimming(QC_ready_reads)
	quality_trimming.out.view()

	remove_polyG(quality_trimming.out)
	remove_polyG.out.view()

	reads_quality_control_post(remove_polyG.out)

	//***************************************************************
	//Reads alignment

	reads_mapping(prepare_genome_bwa.out,remove_polyG.out)

	//***************************************************************
	//Post processing
	
	post_alignment_processing(reads_mapping.out)

	
}

