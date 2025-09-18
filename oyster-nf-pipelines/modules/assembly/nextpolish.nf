process nextpolish {
    tag "${contigs}"
    container 'quay.io/biocontainers/nextpolish:1.4.1--py39hdfd78af_0'

    cpus 16
    memory '64 GB'
    time '24h'

    input:
    path contigs
    path illumina_reads

    output:
    path "polished.fa"

    script:
    """
    nextPolish run -t ${task.cpus} -g ${contigs} -1 ${illumina_reads}
    mv genome.nextpolish.fa polished.fa
    """
}