params.genome_dir="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp"
params.genome="primary_dedup_chr_masked_hp_sealed.fa"
params.cohort_name="20cv_df4_hoh"

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


process gatk_variant_calling{
        container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

        cpus 2
        memory '19 GB'

        input:
        tuple val(replicateID),path(bamfile),path(bamindex)
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        tuple val(replicateID),path("${replicateID}.g.vcf")

        script:
        """
        gatk --java-options "-Xmx18G -Djava.io.tmpdir=/home/zl436/palmer_scratch/simu/tool/nextflow" HaplotypeCaller -R ${ref_fasta} -I ${bamfile} -O ${replicateID}.g.vcf -ERC GVCF
        """

}

process gatk_joint_genotyping{
        container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

        cpus 2
        memory '18 GB'

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

workflow{

        //***************************************************************
        //Generating genome indices

        prepare_genome_samtools("${params.genome_dir}/${params.genome}")

        prepare_genome_picard("${params.genome_dir}/${params.genome}")

        //prepare_genome_bwa("${params.genome_dir}/${params.genome}")

        extract_chromosome_info("${params.genome_dir}/${params.genome}")

        intvl_ch=extract_chromosome_info.out.splitText().map{it->it.trim()}
        intvl_ch.view()

	
        //***************************************************************
        //Variant calling

	Channel.fromPath("/home/zl436/palmer_scratch/bam/2nd_draft_final/*.bam").map { bam_file ->
            def sample_id = bam_file.getBaseName().replaceFirst(/\.bam$/, '')
            def bai_file = file("${bam_file}.bai")
            if (!bai_file.exists()) {
                bai_file = file("${bam_file.getParent()}/${sample_id}.bai")
            }
            tuple(sample_id, bam_file, bai_file)
        }.filter { sample_id, bam, bai -> bai.exists() }.set { bam_ch }

	bam_ch.view()

        gatk_variant_calling(bam_ch,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out)
	gatk_variant_calling.out.view()

        //Generate sample map

        sample_map = gatk_variant_calling.out.collectFile(){id,gvcf->["${params.cohort_name}_map.tsv","${id}\t${gvcf}\n"]}.first()
	sample_map.view()

        //***************************************************************
        //Joint variant calling

        gatk_joint_genotyping(params.cohort_name,sample_map,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out,intvl_ch)
        gatk_joint_genotyping.out.view()


}



