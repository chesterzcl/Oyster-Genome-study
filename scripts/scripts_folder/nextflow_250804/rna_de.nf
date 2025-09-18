#!/usr/bin/env nextflow

// ===============================
// Parameters
// ===============================
params.reads        = "/home/zl436/palmer_scratch/fastq/rna/immune/*_{1,2}.fastq"
params.genome_fasta = "/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
params.gtf_file     = "/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc_hp/gene_pred/braker.gtf"
params.outdir       = "/home/zl436/palmer_scratch/rna_aligned"
params.threads      = 8

// ===============================
// Channel for FASTQ files
// ===============================
Channel
    .fromFilePairs( "${params.reads}", checkIfExists: true )
    .set { read_pairs }

// Create a value channel for genome fasta and gtf
genome_fasta_ch = Channel.value(file(params.genome_fasta))
gtf_file_ch     = Channel.value(file(params.gtf_file))

// ===============================
// Step 0: Build STAR index (runs once)
// ===============================
process build_star_index {
    cpus 8
    memory '40 GB'
    time '12h'

    input:
    path genome_fasta
    path gtf_file

    output:
    path "STARindex"

    script:
    """
    module load STAR
    mkdir -p STARindex
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir STARindex \
         --genomeFastaFiles ${genome_fasta} \
         --sjdbGTFfile ${gtf_file} \
         --sjdbOverhang 149
    """
}

// ===============================
// Step 1: Quality trimming with fastp
// ===============================
process fastp_trim {
    cpus 8
    memory '40 GB'
    time '6h'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_R1.trim.fastq.gz"), path("${sample}_R2.trim.fastq.gz")

    script:
    """
    module load fastp
    fastp \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${sample}_R1.trim.fastq.gz -O ${sample}_R2.trim.fastq.gz \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --thread ${task.cpus} \
        --html ${sample}_fastp.html \
        --json ${sample}_fastp.json
    """
}

// ===============================
// Step 2: STAR alignment
// ===============================
process star_align {
    cpus 8
    memory '40 GB'
    time '12h'

    input:
    tuple val(sample), path(trim_r1), path(trim_r2)
    path star_index

    output:
    tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam")

    script:
    """
    module load STAR
    mkdir -p ${sample}_STAR
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${star_index} \
         --readFilesIn ${trim_r1} ${trim_r2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 30000000000 \
         --outFileNamePrefix ${sample}_STAR/
    mv ${sample}_STAR/Aligned.sortedByCoord.out.bam ${sample}.Aligned.sortedByCoord.out.bam
    """
}

// ===============================
// Step 3: FeatureCounts quantification
// ===============================
process featureCounts {
    cpus 8
    memory '40 GB'
    time '6h'

    input:
    tuple val(sample), path(bam)

    output:
    path "${sample}_counts.txt"

    script:
    """
    module load Subread
    featureCounts -T ${task.cpus} \
                  -a ${params.gtf_file} \
                  -t exon -g gene_id \
                  -p -B -C \
                  -o ${sample}_counts.txt \
                  ${bam}
    """

}

// ===============================
// Workflow Definition
// ===============================

workflow {

    // Step 0: Generate STAR index
    star_index_ch = build_star_index(genome_fasta_ch, gtf_file_ch)

    // Step 1: Trim reads
    trimmed_ch = read_pairs | fastp_trim

    // Step 2: Align reads (pass both channels directly)
    aligned_ch = star_align(trimmed_ch, star_index_ch)

    // Step 3: Quantification
    aligned_ch | featureCounts
}


