params.input_fastq ="/home/zl436/palmer_scratch/fastq/rna/sra/*_{1,2}.fastq"

params.dir="/home/zl436/palmer_scratch/oyster_genome"
params.reads_wgs ="${params.dir}/WGS/*_{1,2}.fastq.gz"
params.genome_fasta ="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
params.ref_fasta="/home/zl436/palmer_scratch/ref/CV.fa"
params.repeat_lib="/home/zl436/palmer_scratch/oyster_genome/chr_level_asmbl/custom_repeats_dedup.fa"
params.gtf = "/home/zl436/palmer_scratch/oyster_genome/Contigs_scfd_salsa/annotation.gtf"
params.outdir="/home/zl436/palmer_scratch/oyster_genome/Contigs_scfd_salsa/rna_bams"
params.threads =4

process QC_RAW {
	container 'quay.io/biocontainers/fastqc:0.11.9--0'

        cpus 4
        memory '50GB'

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample),path("${sample}_*_fastqc.zip")

	script:
	"""
	fastqc -t ${params.threads} ${reads}
	"""
}


process TRIM_FASTQ {
	container 'quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0'

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}_trimmed*")

	script:
	"""
	fastp -i ${sample}_1.fastq -I ${sample}_2.fastq -o ${sample}_trimmed_1.fastq.gz -O ${sample}_trimmed_2.fastq.gz --cut_front --cut_tail --detect_adapter_for_pe  --cut_mean_quality 20 -l 50
	"""
}

process QC_TRIMMED {
	container 'quay.io/biocontainers/fastqc:0.11.9--0'

        cpus 4
        memory '50GB'

	input:
	tuple val(sample), path(trimmed_reads)

	output:
	tuple val(sample),path("${sample}_*_fastqc.zip")

	script:
	"""
	fastqc -t ${params.threads} ${trimmed_reads}
	"""
}

process evaluate_contigs{
        container 'quay.io/biocontainers/quast:5.3.0--py311pl5321hc84137b_1'

        cpus 8

        input:
        tuple val(assembly_name),path(fa_files)

        output:
        tuple val(assembly_name),path("${assembly_name}_quast_report")

        script:
        """
        quast -t ${task.cpus} -o ${assembly_name}_quast_report ${fa_files}
        """

}

process evaluate_genes{
        container 'quay.io/biocontainers/busco:5.8.2--pyhdfd78af_0'

        cpus 8

        input:
        tuple val(assembly_name),path(fa_files)

        output:
        tuple val(assembly_name),path("${assembly_name}_busco_report")

        script:
        """
        busco -c ${task.cpus} -m genome -i ${fa_files} -o ${assembly_name}_busco_report -l mollusca_odb10
        """

}

process get_kmer_content{
        container 'quay.io/biocontainers/meryl:1.4.1--h9948957_2'

        cpus 8

        input:
        tuple val(sample_id),path(reads1),path(reads2)

        output:
        path("${sample_id}_kmer_db_filt_1.meryl")

        script:
        """
        meryl k=21 count output ${sample_id}_kmer_db.meryl ${reads1} ${reads2}
	meryl greater-than 1 ${sample_id}_kmer_db.meryl ${sample_id}_kmer_db_filt_1.meryl
        """

}

process evaluate_kmers{

        container 'danylmb/merqury:1.4.1-build2'

        cpus 8

        input:
        tuple val(assembly_name),path(fa_files),path(kmer_db)

        output:
        tuple val(assembly_name),path("${assembly_name}_merqury_output*")

        script:
        """
        merqury.sh ${kmer_db} ${fa_files} ${assembly_name}_merqury_output
        """
}

process repeat_modeling{
	container 'quay.io/biocontainers/repeatmodeler:1.0.11--pl5.22.0_0'

	cpus 8

	input:
	tuple val(assembly_name),path(fa_file)

	output:
	tuple val(assembly_name),path(fa_file),path("*.classified")
	
	script:
	"""
	BuildDatabase -name oyster_db ${fa_file}
	RepeatModeler -database oyster_db -pa 8
	"""
}

process repeat_masking{
	container 'quay.io/biocontainers/repeatmasker:4.1.8--pl5321hdfd78af_0'

	cpus 8

	input:
	tuple val(assembly_name),path(fa_file),path(reapeat_db)

	output:
	tuple val(assembly_name),path("*.fa.masked")

	script:
	"""
	RepeatMasker -pa 8 -lib ${repeat_db} -xsmall ${fa_file}
	
	"""

}

process repeat_masking_pri{
        container 'quay.io/biocontainers/repeatmasker:4.1.8--pl5321hdfd78af_0'

        cpus 16

        input:
        tuple val(assembly_name),path(fa_file)

        output:
        tuple val(assembly_name),path("*.fa.masked")

        script:
        """
        RepeatMasker -pa 16 -lib ${params.repeat_lib} -xsmall ${fa_file}

        """
}

process STAR_INDEX {
	container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_5'

        cpus 4
        memory '50GB'
	
	input:
	tuple val(assembly_name),path(fasta)

	output:
	path genome_index

	script:
	"""
	STAR --runThreadN ${params.threads} \
	--runMode genomeGenerate \
	--genomeDir genome_index \
	--genomeFastaFiles ${fasta} 
	"""
}


process STAR_ALIGN {
	container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_5'

	cpus 4
	memory '50GB'

	input:
	tuple val(sample),path(trimmed_reads),path(genome_index)

	output:
	tuple val(sample),path("${sample}_Aligned.sortedByCoord.out.bam")

	script:
	"""
	STAR --runThreadN ${task.cpus} \
	--genomeDir genome_index \
	--readFilesIn ${trimmed_reads} \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${sample}_ \
	--limitBAMsortRAM 20000000000

	"""
}

process postprocessing{
	container 'quay.io/biocontainers/samtools:1.6--h5fe306e_11'

	input:
	tuple val(sample),path(bam_file)

	output:
	tuple val(sample),path(bam_file),path("*.bam.bai")

	script:
	"""
	samtools index ${bam_file}
	samtools flagstat ${bam_file}
	"""
}

process merge_bams{
        container 'quay.io/biocontainers/samtools:1.6--h5fe306e_11'

	cpus 8

	input:
	path bams
	path indices

	output:
	tuple path("rna_merged.bam"),path("rna_merged.bam.bai")

	script:
	"""
	samtools merge -@ 8 rna_merged.bam ${bams}
	samtools index rna_merged.bam
	"""

}

process predict_genes{
	container 'teambraker/braker3'

	cpus 16

	input:
	tuple val(assembly_name),path(masked_genome)
	tuple path(merged_bam),path(bam_index)

	output:
	path "braker_output/*"

	script:
	"""
	braker.pl --genome=${masked_genome} --bam=${merged_bam} --softmasking --threads=16 --species=oyster --workingdir=braker_output
	"""

}

process transcript_assembly{
	container 'quay.io/biocontainers/stringtie:3.0.0--h29c0135_0'

	input:
	tuple val(sample),path(bam_file),path(bam_idx)

	output:
	tuple val(sample),path("*.gtf")

	script:
	"""
	stringtie ${bam_file} -o ${sample}.gtf
	"""
}

process merge_transcript{
	container 'quay.io/biocontainers/stringtie:3.0.0--h29c0135_0'

	input:
	path gtfs
	
	output:
	path("merged_transcripts.gtf")

	script:
	"""
	stringtie --merge -o merged_transcripts.gtf ${gtfs.join(' ')}
	"""
}

process compare_genome{
	container 'quay.io/biocontainers/mummer4:4.0.0--pl5321h9948957_0'

	input:
	tuple val(ref_name),path(genome_ref)
	tuple val(scfd_name),path(genome_scfd)

	output:
	tuple path("${ref_name}_${scfd_name}.coords"),path("${ref_name}_${scfd_name}.delta")

	script:
	"""
	nucmer --prefix=${ref_name}_${scfd_name} ${genome_ref} ${genome_scfd}
	mummerplot --png --layout --prefix=${ref_name}_${scfd_name} ${ref_name}_${scfd_name}.delta
	show-coords -rcl ${ref_name}_${scfd_name}.delta>${ref_name}_${scfd_name}.coords
	"""
}


process gene_prediction {
	container 'quay.io/biocontainers/augustus:3.5.0--pl5321heb9362c_5'

	input:
	path merged_transcripts
	tuple val(genome_name),path(genome_fasta)

	output:
	path "predicted_genes.gff"

	script:
	"""
	augustus --genome=${genome_fasta} --strand=both --softmasking=on --exonnames=on --gtf=1 --output=predicted_genes.gff
	"""
}

workflow {
	samples = Channel.fromFilePairs(params.input_fastq)
	samples.view()	

	Channel.fromFilePairs(params.reads_wgs,flat:true).set{wgs_ch}
	//kmer_ch=wgs_ch|get_kmer_content
        //kmer_ch.view()

	genome=Channel.fromPath(file(params.genome_fasta)).map{file->tuple(file.baseName,file)}
	genome.view()

        genome|evaluate_contigs
        evaluate_contigs.out.view()

        genome|evaluate_genes
        evaluate_genes.out.view()

        //comb_kmer_ch=genome.combine(kmer_ch)
        //comb_kmer_ch.view()

        //evaluate_kmers(comb_kmer_ch)
        //evaluate_kmers.out.view()

	samples|QC_RAW
	samples|TRIM_FASTQ|QC_TRIMMED
    
	//masked_genome=genome|repeat_modeling|repeat_masking
	//masked_genome=genome|repeat_masking_pri
	//masked_genome.view()

	genome_idx_ch=genome|STAR_INDEX
	trimmed_samples_ch=TRIM_FASTQ.out

    	combined_samples_ch = trimmed_samples_ch.combine(genome_idx_ch).map {[it[0],it[1],it[2]]}
	combined_samples_ch.view()

	aligned_rna_ch=combined_samples_ch|STAR_ALIGN
	aligned_rna_ch.view()

	processed_bam=aligned_rna_ch|postprocessing
	processed_bam.view()

	bam_modified_ch=processed_bam.map{it[1]}.toList()
	bam_modified_ch.view()	

	bam_index_ch=processed_bam.map{it[2]}.toList()
	bam_index_ch.view()

	merged_bam=merge_bams(bam_modified_ch,bam_index_ch)
	merged_bam.view()

	predicted_ch=predict_genes(genome,merged_bam)
	predicted_ch.view()

	//gtf_ch=aligned_rna_ch|postprocessing|transcript_assembly
	//gtf_ch.view()
	//gtf_modified_ch=gtf_ch.map{it[1]}.toList()
	//gtf_modified_ch.view()

	//merged_gtf_ch=gtf_modified_ch|merge_transcript
	//merged_gtf_ch.view()

	//ref=Channel.fromPath(file(params.ref_fasta)).map{file->tuple(file.baseName,file)}
	//compare_genome(ref,genome)
	//compare_genome.out.view()	

}














