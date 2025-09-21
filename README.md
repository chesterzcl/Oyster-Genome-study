# Multi-omic Investigation of the Eastern Oyster Genome


Chromosome-level assembly and multi-omics analysis of the eastern oyster (*Crassostrea virginica*).
This repository contains code/pipeline of the genome assembly, annotation, and comparative/evolutionary analyses.


*Schematic overview of the study*

![Assembly workflow schematic](docs/figures/Schematic.png)


## Contents
- **Assembly** — PacBio HiFi + Illumina + Hi-C integration (hifiasm, purge_dups, 3D-DNA, Juicebox curation)
- **Annotation** — BRAKER gene prediction, repeat annotation, regulatory landscapes
- **Custom scripts** — integrative analyses (SNP↔ATAC overlap, Circos plots, synteny parsing)
- **Evolutionary analysis** — synteny alignments, chromosome rearrangements, gene family evolution