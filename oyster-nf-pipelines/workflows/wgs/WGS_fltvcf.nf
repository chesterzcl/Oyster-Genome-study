params.genome_dir="/vast/palmer/scratch/hoh/zl436/oyster_genome/Contigs_dedup_scfd_3ddna_mc"
params.genome="primary_dedup_chr.fa"
params.vcf_dir="/home/zl436/palmer_scratch/vcf/oyster/snp/2nd_draft_new"
params.vcf="${params.vcf_dir}/20cv_df3_hoh.vcf"
params.vcf_idx="${params.vcf_dir}/20cv_df3_hoh.vcf.idx"

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


process gatk_vcf_filtering{
        container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

        cpus 2
        memory '11 GB'

        input:
	path vcf
	path vcf_idx
        path ref_fasta
        path ref_index
        path ref_dict

        output:
	tuple path("*_filtered.vcf"),path("*_filtered.vcf.idx")

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


workflow{
        //***************************************************************
        //Generating genome indices
        prepare_genome_samtools("${params.genome_dir}/${params.genome}")
        prepare_genome_picard("${params.genome_dir}/${params.genome}")

	//***************************************************************
	//Read files
	vcf_ch=Channel.fromPath(params.vcf,checkIfExists: true)
	vcf_idx_ch=Channel.fromPath(params.vcf_idx,checkIfExists: true)
	
	//***************************************************************
	//SNP_filtereing
	gatk_vcf_filtering(vcf_ch,vcf_idx_ch,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out)
	gatk_vcf_filtering.out.view()


}


















