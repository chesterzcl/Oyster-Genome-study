process busco {
    tag "${scaffolds}"
    container 'ezlabgva/busco:v5.4.4_cv1'

    cpus 8
    memory '32 GB'

    input:
    path scaffolds

    output:
    path "busco_results/"

    script:
    """
    busco -i ${scaffolds} -l mollusca_odb10 -o busco_results -m genome -c ${task.cpus}
    """
}