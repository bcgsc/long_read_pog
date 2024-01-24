# The Long-read POG Cohort
[![render-rmarkdown](https://github.com/bcgsc/nanopore_pog/actions/workflows/renderRMarkdown.yaml/badge.svg)](https://github.com/bcgsc/nanopore_pog/actions/workflows/renderRMarkdown.yaml)

The Long-read POG Cohort is a set of whole-genome sequencing of diverse advanced cancer samples from participants in the [Personalized Oncogenomics Programme at Canada's Michael Smith Genome Sciences Centre](https://www.bcgsc.ca/personalized-oncogenomics-program), sequenced using long-read sequencing modalities. The first release includes 189 tumours from 181 individuals sequenced on an Oxford Nanopore Technologies PromethION, using the R9.4.1 chemistry. 41 of these include matching normal samples from the same individuals. 

## Raw sequence data and germline variants

## Metadata and clinical information


## Figure generation code
All code to generate all figures in the main manuscript can be found under the `FigureX/` subdirectories. All code uses anonymised summary data made available at [https://www.bcgsc.ca/downloads/nanopore_pog](https://www.bcgsc.ca/downloads/nanopore_pog/), and should be possible to run directly.


## Analysis code
Code for the more computationally-intensive analyses that created the summary data is provided under the `Analysis/` subdirectories. This will generally require some configuration (and access to HPC resources) to get working. This will also require downloading the raw sequence data from EGA. The primary purpose of providing this code is to serve as a record of tool versions and parameters used.
