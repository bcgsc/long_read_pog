# The Long-read POG Cohort
[![render-rmarkdown](https://github.com/bcgsc/nanopore_pog/actions/workflows/renderRMarkdown.yaml/badge.svg?branch=main)](https://github.com/bcgsc/nanopore_pog/actions/workflows/renderRMarkdown.yaml)
[![](https://img.shields.io/docker/image-size/bcgsc/long-pog-rmd?label=Rmd%20Docker%20image)](https://hub.docker.com/repository/docker/bcgsc/long-pog-rmd)
[![](https://img.shields.io/docker/image-size/bcgsc/long-pog-jupyter?label=Jupyter%20Docker%20image)](https://hub.docker.com/repository/docker/bcgsc/long-pog-jupyter)

The Long-read POG Cohort is a set of whole-genome sequencing of diverse advanced cancer samples from participants in the [Personalized Oncogenomics Programme at Canada's Michael Smith Genome Sciences Centre](https://www.bcgsc.ca/personalized-oncogenomics-program), sequenced using long-read sequencing modalities. The first release includes 189 tumours from 181 individuals sequenced on an Oxford Nanopore Technologies PromethION, using the R9.4.1 chemistry. 41 of these include matching normal samples from the same individuals. 


## Metadata and clinical information
Sample information, high-level clinical information, as well as final data used for figures is available at [https://www.bcgsc.ca/downloads/nanopore_pog](https://www.bcgsc.ca/downloads/nanopore_pog/).

## Raw sequence data and germline variants
As this data contains sensitive information, it is securely stored at the [European Genome-Pheneom Archive (EGA) under accession EGAS00001001159](https://ega-archive.org/studies/EGAS00001001159). That study accession also includes Illumina data for these and other POG cases. Per-sample accessions are listed in the [main sample sheet](https://www.bcgsc.ca/downloads/nanopore_pog/supplementary_tables/Supplementary_Table_1_samples.tsv). 

This data is broadly consented for use in research, with some conditions, such as proof of IRB approval. Access requests via EGA will go to our data access committee, who will ensure that the request meets ethical requirements. 


## Figure generation code
All code to generate all figures in the manuscript can be found under `Figures/` in this repository. Some exceptions exist for figures that were generated using tools such as IGV, and some figures were lightly post-processed in Inkscape or Illustrator (without altering the data representation).

This code will load data from bcgsc.ca. Docker containers are provided, enabling the figures to be reproduced in any environment that can run Docker or Singularity. 

## Analysis code
Code for the more computationally-intensive analyses that created the summary data is provided under the `Analysis/` subdirectories. This will generally require some configuration (and access to HPC resources) to get working. This will also require downloading the raw sequence data from EGA. The primary purpose of providing this code is to serve as a record of tool versions and parameters used.
