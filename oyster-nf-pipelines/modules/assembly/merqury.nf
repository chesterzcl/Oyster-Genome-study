process merqury {
    tag "${scaffolds}"
    container 'quay.io/biocontainers/merqury:1.3--hdfd78af_0'

    cpus 8
    memory '32 GB'

    input:
    path scaffolds
    path hifi_reads

    output:
    path "merqury_out/"

    script:
    """
    merqury.sh ${hifi_reads} ${scaffolds} merqury_out
    """
}