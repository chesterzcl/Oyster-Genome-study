process salsa {
    tag "${polished}"
    container 'quay.io/biocontainers/salsa:2.2--py36h516909a_1'

    cpus 16
    memory '64 GB'
    time '24h'

    input:
    path polished
    path hic_reads

    output:
    path "scaffolds.fa"

    script:
    """
    run_pipeline.py -a ${polished} -l assembly.fai -b hic.bam -o out_salsa
    cp out_salsa/scaffolds_FINAL.fasta scaffolds.fa
    """
}