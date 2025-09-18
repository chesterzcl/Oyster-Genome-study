#!/bin/sh
set -xve
hostname
cd /vast/palmer/scratch/hoh/zl436/simu/tool/nextflow_ga/work/ff/d41971e123aff44462a6d707cf9554/01_rundir/01.hifi_polish/04.merge.bam.sh.work/merge_bam1
time /opt/NextPolish/bin/samtools merge -f -b /vast/palmer/scratch/hoh/zl436/simu/tool/nextflow_ga/work/ff/d41971e123aff44462a6d707cf9554/./01_rundir/01.hifi_polish/hifi.sort.bam.list --threads 2 /vast/palmer/scratch/hoh/zl436/simu/tool/nextflow_ga/work/ff/d41971e123aff44462a6d707cf9554/./01_rundir/01.hifi_polish/hifi.sort.bam
time /opt/NextPolish/bin/samtools index -@ 5 /vast/palmer/scratch/hoh/zl436/simu/tool/nextflow_ga/work/ff/d41971e123aff44462a6d707cf9554/./01_rundir/01.hifi_polish/hifi.sort.bam
touch /vast/palmer/scratch/hoh/zl436/simu/tool/nextflow_ga/work/ff/d41971e123aff44462a6d707cf9554/01_rundir/01.hifi_polish/04.merge.bam.sh.work/merge_bam1/nextPolish.sh.done

