params.genome_dir="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp"
params.genome="primary_dedup_chr_masked_hp_sealed.fa"
params.cohort_name="20cv_final_hoh"
params.data_dir="/home/zl436/ycga_work/oyster_genome/WGS"
params.reads="${params.data_dir}/*_{1,2}.fastq.gz"


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

process extract_chromosome_info{
        input:
        path genome

        output:
        path "chromosome_list.txt"

        script:
        """
        grep '^>' ${genome}|cut -d " " -f 1|cut -c 2- >"chromosome_list.txt"
        """
}

process reads_quality_control{

        container 'quay.io/biocontainers/fastqc:0.11.9--0'

        cpus 4
	memory '40GB'

        input:
        tuple val(replicateID),path(reads)

        output:
        tuple val(replicateID),path("*_1_fastqc.zip"),path("*_2_fastqc.zip")

        script:
        """
        fastqc -t ${task.cpus} ${reads}
        """
}

process reads_quality_control_cp1{

        container 'quay.io/biocontainers/fastqc:0.11.9--0'

        cpus 4

        input:
        tuple val(replicateID),path(reads)

        output:
        tuple val(replicateID),path("*_1_trimmed_grmv_fastqc.zip"),path("*_2_trimmed_grmv_fastqc.zip")

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
        trimming_params="ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:15 TRAILING:5 SLIDINGWINDOW:10:30 MINLEN:50"
        """
}

process quality_trimming{

        container 'quay.io/biocontainers/trimmomatic:0.36--3'

        cpus 4

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
        fastp -i ${replicateID}_1_trimmed.fastq.gz -I ${replicateID}_2_trimmed.fastq.gz -o ${replicateID}_1_trimmed_grmv.fastq.gz -O ${replicateID}_2_trimmed_grmv.fastq.gz --trim_poly_g -l 50 --thread ${task.cpus}
        """
}

process reads_mapping{
        container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

        cpus 4

        input:
        path genome_index
        tuple val(replicateID),path(trimmed_reads)

        output:
        tuple val(replicateID),path("*.sam")

        script:
        """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${replicateID}\\tSM:${replicateID}.SM\\tLB:${replicateID}.LB' -o ${replicateID}.sam ${genome_index[0].baseName} ${trimmed_reads}
        rm ${trimmed_reads}
        """
}

process post_alignment_processing{
        container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

        cpus 4

        input:
        tuple val(replicateID),path(samfile)

        output:
        tuple val(replicateID),path("${replicateID}_st_dr_filtered.bam"),path("${replicateID}_st_dr_filtered.bam.bai")

        script:
        """
        samtools fixmate -@ ${task.cpus} -O BAM ${replicateID}.sam ${replicateID}.bam
        rm ${replicateID}.sam
        samtools sort -@ ${task.cpus} -O BAM -o ${replicateID}_sorted.bam ${replicateID}.bam
        rm ${replicateID}.bam
        samtools rmdup ${replicateID}_sorted.bam ${replicateID}_st_dr.bam
        rm ${replicateID}_sorted.bam
        samtools view -@ ${task.cpus} -O BAM -q 30 -o ${replicateID}_st_dr_filtered.bam ${replicateID}_st_dr.bam
        rm ${replicateID}_st_dr.bam
        samtools index ${replicateID}_st_dr_filtered.bam
        """
}

process gatk_variant_calling{
        container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

        cpus 2
        memory '20GB'

        input:
        tuple val(replicateID),path(bamfile),path(bamindex)
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        tuple val(replicateID),path("${replicateID}.g.vcf")

        script:
        """
        gatk HaplotypeCaller -R ${ref_fasta} -I ${bamfile} -O ${replicateID}.g.vcf -ERC GVCF
        """

}

process gatk_joint_genotyping{
        container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

        cpus 2
        memory '20GB'

        input:
        val cohort_name
        path sample_map
        path ref_fasta
        path ref_index
        path ref_dict
        val interval_name

        output:
        tuple val(interval_name),path("${cohort_name}_${interval_name}.vcf"),path("${cohort_name}_${interval_name}.vcf.idx")

        script:
        """
        gatk GenomicsDBImport --sample-name-map ${sample_map} --genomicsdb-workspace-path ${cohort_name}_${interval_name}_db -L ${interval_name}

        gatk GenotypeGVCFs -R ${ref_fasta} -V gendb://${cohort_name}_${interval_name}_db -O ${cohort_name}_${interval_name}.vcf -L ${interval_name}

        """
}

process gatk_vcf_filtering{
        container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

        cpus 2
        memory '11 GB'

        input:
        tuple val(interval_name),path(vcf),path(vcf_idx)
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        tuple val(interval_name),path("*_filtered.vcf"),path("*_filtered.vcf.idx")

        script:
        """
        gatk SelectVariants --restrict-alleles-to BIALLELIC --select-type-to-include SNP -R ${ref_fasta} -V ${vcf} -O ${vcf.baseName}_biSNP.vcf
        gatk IndexFeatureFile -I ${vcf.baseName}_biSNP.vcf

        gatk VariantFiltration -V ${vcf.baseName}_biSNP.vcf -O ${vcf.baseName}_biSNP_marked.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
        gatk IndexFeatureFile -I ${vcf.baseName}_biSNP_marked.vcf

        gatk SelectVariants --restrict-alleles-to BIALLELIC --exclude-filtered --select-type-to-include SNP -R ${ref_fasta} -V ${vcf.baseName}_biSNP_marked.vcf -O ${vcf.baseName}_biSNP_filtered.vcf
        gatk IndexFeatureFile -I ${vcf.baseName}_biSNP_filtered.vcf

        """
}

process merge_vcf {
	container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

    	input:
    	tuple val(chr_name), path(vcf_file), path(vcf_index)
	val cohort_name

	output:
    	path "${cohort_name}.vcf"

    	script:
    	"""
	gatk MergeVcfs \
	-I ${vcf_file.join(' -I ')} \
	-O ${cohort_name}.vcf

	gatk IndexFeatureFile -I ${cohort_name}.vcf
	"""
}

workflow{

        //***************************************************************
        //Generating genome indices

        prepare_genome_samtools("${params.genome_dir}/${params.genome}")

        prepare_genome_picard("${params.genome_dir}/${params.genome}")

        prepare_genome_bwa("${params.genome_dir}/${params.genome}")

        extract_chromosome_info("${params.genome_dir}/${params.genome}")

        //***************************************************************
        //Prepare input data

        reads_ch=Channel.fromFilePairs(params.reads,checkIfExists: true)
        reads_ch.view()
        intvl_ch=extract_chromosome_info.out.splitText().map{it->it.trim()}
        intvl_ch.view()

        //***************************************************************
        //Preprocessing

        reads_quality_control(reads_ch)
	reads_quality_control.out.view()

        get_trimming_parameters(reads_quality_control.out)
        get_trimming_parameters.out.view()

        reads_ch=Channel.fromFilePairs(params.reads,checkIfExists: true)
        reads_ch.join(get_trimming_parameters.out).set{QC_ready_reads}
        QC_ready_reads.view()

        quality_trimming(QC_ready_reads)
        quality_trimming.out.view()

        remove_polyG(quality_trimming.out)
	trimmed_reads_ch=remove_polyG.out
        trimmed_reads_ch.view()
	
	trimmed_qc_ch=trimmed_reads_ch|reads_quality_control_cp1
	trimmed_qc_ch.view()
        //***************************************************************
        //Reads alignment
	
        reads_mapping(prepare_genome_bwa.out,trimmed_reads_ch)
	reads_mapping.out.view()

        //***************************************************************
        //Post processing

        post_alignment_processing(reads_mapping.out)
	post_alignment_processing.out.view()	

        //***************************************************************
        //Variant calling

        //gatk_variant_calling(post_alignment_processing.out,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out)
	//gatk_variant_calling.out.view()

        //Generate sample map

        //sample_map = gatk_variant_calling.out.collectFile(){id,gvcf->["${params.cohort_name}_map.tsv","${id}\t${gvcf}\n"]}.first()
	//sample_map.view()

        //***************************************************************
        //Joint variant calling

        //raw_vcfs=gatk_joint_genotyping(params.cohort_name,sample_map,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out,intvl_ch)
        //raw_vcfs.view()

	//Filter vcfs

	//filtered_vcfs=gatk_vcf_filtering(raw_vcfs,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out)
	//filtered_vcfs.view()	

	//merged_vcf=merge_vcf(filtered_vcfs,params.cohort_name)
	//merged_vcf.view()
}



