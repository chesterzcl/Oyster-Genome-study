include { hifiasm }   from '../modules/assembly/hifiasm'
include { nextpolish } from '../modules/assembly/nextpolish'
include { salsa }     from '../modules/assembly/salsa'
include { busco }     from '../modules/assembly/busco'
include { merqury }   from '../modules/assembly/merqury'

workflow assembly {
    take:
    hifi_reads
    illumina_reads
    hic_reads

    main:
    contigs   = hifiasm(hifi_reads)
    polished  = nextpolish(contigs, illumina_reads)
    scaffolds = salsa(polished, hic_reads)
    qc_busco  = busco(scaffolds)
    qc_merqury = merqury(scaffolds, hifi_reads)

    emit:
    scaffolds
    qc_busco
    qc_merqury
}