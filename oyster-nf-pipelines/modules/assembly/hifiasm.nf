process hifiasm {
    tag "${hifi_reads}"
    container 'quay.io/biocontainers/hifiasm:0.24.0--h5ca1c30_0'

    cpus 16
    memory '64 GB'
    time '24h'

    input:
    path hifi_reads

    output:
    path "contigs.gfa", emit: gfa
    path "contigs.fa",  emit: fasta

    script:
    """
    hifiasm -o contigs -t ${task.cpus} ${hifi_reads}
    # convert to fasta
    awk '/^S/{print ">"\$2"\\n"\$3}' contigs.gfa > contigs.fa
    """
}