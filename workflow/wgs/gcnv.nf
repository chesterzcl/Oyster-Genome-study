params.bams = "/home/zl436/palmer_scratch/bam/2nd_draft/*.bam"
params.ref = "/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/primary_dedup_chr.fa"
params.bin_length = 500
params.repeat="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/repeat_analysis/cnv_masked_repeats.bed"

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

process prepare_genome_samtools {
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

process PreprocessIntervals {
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_0'

    cpus 2
    memory '20 GB'

    input:
    path fasta
    path dict
    path fai

    output:
    path "bins.interval_list"

    script:
    """
    gatk PreprocessIntervals \
	-R ${fasta} \
	--bin-length ${params.bin_length} \
	--padding 0 \
	--interval-merging-rule OVERLAPPING_ONLY \
	-O bins.interval_list
    """
}

process PreFilterIntervals {
    container 'quay.io/biocontainers/bedtools:2.27.1--h077b44d_9'

    input:
    path intervals_list         // From PreprocessIntervals (.interval_list)
    path repeat_bed             // From params.repeat or a previous process

    output:
    path "filtered_intervals.interval_list"  // Ready for GATK downstream

    script:
    """
    awk 'BEGIN{OFS="\\t"} !/^@/ {print \$1, \$2-1, \$3}' ${intervals_list} > pre_filter_intervals.bed

    bedtools subtract -a pre_filter_intervals.bed -b ${repeat_bed} > filtered.bed

    grep '^@' ${intervals_list} > filtered_intervals.interval_list
    awk 'BEGIN{OFS="\\t"} {print \$1,\$2+1,\$3,"+","N"}' filtered.bed >> filtered_intervals.interval_list
    """
}

process RebinMaskedIntervals {
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_0'

    cpus 2
    memory '16 GB'
    
    input:
    path fasta
    path dict
    path fai
    path masked_intervals_list

    output:
    path "rebinned_intervals.interval_list"

    script:
    """
    gatk PreprocessIntervals \\
        -R ${fasta} \\
        -L ${masked_intervals_list} \\
        --bin-length ${params.bin_length} \\
        --padding 0 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O rebinned_intervals.interval_list
    """
}


process AnnotateIntervals {
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_0'
    cpus 2
    memory '10 GB'

    input:
    path interval_list
    path fasta
    path fai
    path dict

    output:
    path "annotated_intervals.tsv"

    script:
    """
    gatk AnnotateIntervals \
        -R ${fasta} \
        -L ${interval_list} \
        -O annotated_intervals.tsv \
    	--interval-merging-rule OVERLAPPING_ONLY
    """
}

process CollectReadCounts {
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_0'
    cpus 2
    memory '10 GB'

    input:
    tuple path(bam), path(interval)

    output:
    path "${bam.simpleName}.counts.tsv"

    script:
    """
    gatk CollectReadCounts \
        -I ${bam} \
        -L ${interval} \
        --format TSV \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${bam.simpleName}.counts.tsv
    """
}

process FilterIntervals {
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_0'
    cpus 2
    memory '20 GB'

    input:
    path interval_list
    path annotated_intervals
    path coverage_tsvs

    output:
    path "filtered_cohort_intervals.interval_list"

    script:
    """
    gatk FilterIntervals \
        -L ${interval_list} \
        --annotated-intervals ${annotated_intervals} \
        ${coverage_tsvs.collect { "-I '${it}'" }.join(' ')} \
        -imr OVERLAPPING_ONLY \
        -O filtered_cohort_intervals.interval_list
    """
}


process DeterminePloidy {
    container 'broadinstitute/gatk:4.6.2.0'
    cpus 2
    memory '20 GB'

    input:
    path interval_list
    path count_files

    output:
    path "ploidy-calls/"

    script:
    """
    echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3" > contig-ploidy-priors.tsv
    awk '/^@SQ/ {
        for (i = 1; i <= NF; i++) {
            if (\$i ~ /^SN:/) {
                split(\$i, a, ":");
                contig = a[2];
                print contig "\\t0.01\\t0.01\\t0.97\\t0.01"
            }
        }
    }' ${interval_list} >> contig-ploidy-priors.tsv

    gatk DetermineGermlineContigPloidy \
        -L ${interval_list} \
        --interval-merging-rule OVERLAPPING_ONLY \
        ${count_files.collect { "-I ${it}" }.join(' ')} \
	--contig-ploidy-priors contig-ploidy-priors.tsv \
        --output ploidy-calls \
        --output-prefix cohort_ploidy \
        --verbosity INFO
    """
}


process GermlineCNVCaller {
    container 'broadinstitute/gatk:4.6.2.0'
    cpus 4
    memory '40 GB'

    input:
    path interval_list
    path annotated_intervals
    path ploidy_calls
    path count_files

    output:
    tuple path("cnv-output/cohort_model-model/"), path("cnv-output/cohort_model-calls/")

    script:
    """
    gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L ${interval_list} \
        ${count_files.collect { "-I ${it}" }.join(' ')} \
        --contig-ploidy-calls ${ploidy_calls}/cohort_ploidy-calls \
        --annotated-intervals ${annotated_intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output cnv-output \
        --output-prefix cohort_model \
        --verbosity DEBUG
    """
}

process PostprocessCNVCalls {
    container 'broadinstitute/gatk:4.6.2.0'
    cpus 2
    memory '10 GB'

    input:
    tuple val(sample_index), val(sample_name), path(model_dir), path(call_dir), path(ploidy_calls), path(dict_file)

    output:
    path "genotyped-intervals-${sample_name}.vcf.gz"
    path "genotyped-segments-${sample_name}.vcf.gz"

    script:
    """
    gatk PostprocessGermlineCNVCalls \\
        --model-shard-path ${model_dir} \\
        --calls-shard-path ${call_dir} \\
        --contig-ploidy-calls ${ploidy_calls} \\
        --sample-index ${sample_index} \\
        --output-genotyped-intervals genotyped-intervals-${sample_name}.vcf.gz \\
        --output-genotyped-segments genotyped-segments-${sample_name}.vcf.gz \\
        --sequence-dictionary ${dict_file}
    """
}


workflow {
    ref_ch  = Channel.fromPath(params.ref)
    fai_ch  = ref_ch | prepare_genome_picard
    dict_ch = ref_ch | prepare_genome_samtools
    bam_ch  = Channel.fromPath(params.bams)
    repeat_ch = Channel.fromPath(params.repeat)
    repeat_ch.view()

    // Step 1: Initial binning
    intervals_ch = PreprocessIntervals(ref_ch, dict_ch, fai_ch)

    // Step 2: Mask with repeats.bed
    filtered_intervals_ch = PreFilterIntervals(intervals_ch, repeat_ch)
    filtered_intervals_ch.view()

    // Step 3: Re-bin the repeat-masked intervals
    rebinned_ch = RebinMaskedIntervals(ref_ch, dict_ch, fai_ch, filtered_intervals_ch)
    rebinned_ch.view()

    // Step 4: Annotate rebinned intervals
    annotated_ch = AnnotateIntervals(rebinned_ch, ref_ch, fai_ch, dict_ch)
    annotated_ch.view()

    // Step 5: Collect read counts
    collect_input_ch = bam_ch.combine(rebinned_ch)
    CollectReadCounts(collect_input_ch)
    counts_ch = CollectReadCounts.out.collect()

    // Step 6: Filter intervals and determine ploidy
    filtered_ch = FilterIntervals(rebinned_ch, annotated_ch, counts_ch.flatten())
    filtered_ch.view()

    DeterminePloidy(filtered_ch, counts_ch)
    ploidy_ch = DeterminePloidy.out
    ploidy_ch.view()

    // Step 7: CNV calling
    GermlineCNVCaller(filtered_ch, annotated_ch, ploidy_ch, counts_ch)
    GermlineCNVCaller.out.view()

    // Extract model and call directory paths
    cnv_model_dir_ch = GermlineCNVCaller.out.map { it[0] }
    cnv_calls_dir_ch = GermlineCNVCaller.out.map { it[1] }

    cnv_model_dir_ch.view()
    cnv_calls_dir_ch.view()

    // Generate (index, sample_name) from BAM file base names
    sample_index_ch = bam_ch.map { it.getBaseName() }
                            .collect()
                            .flatMap { names -> names.withIndex().collect { i, name -> tuple(i, name) } }

    sample_index_ch.view()

    // Compose input tuples for PostprocessCNVCalls
//    postprocess_input_ch = sample_index_ch
//    .combine(cnv_model_dir_ch)
//    .combine(cnv_calls_dir_ch)
//    .combine(ploidy_ch)
//    .combine(dict_ch)
//   .map { [[[[index, name], model], call], ploidy], dict ->
//        tuple(index, name, model, call, ploidy, dict)}
    
//    postprocess_input_ch.view()


    // Run the postprocessing step
//    PostprocessCNVCalls(postprocess_input_ch)
//    PostprocessCNVCalls.out.view()

}


