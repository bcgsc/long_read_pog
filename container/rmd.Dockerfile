FROM bioconductor/bioconductor_docker:RELEASE_3_17-R-4.3.0

# R packages:
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("data.table", "1.14.8")'
RUN R -e 'remotes::install_version("tidyverse", "2.0.0")'
RUN R -e 'remotes::install_version("MetBrewer", "0.2.0")'
RUN R -e 'remotes::install_version("ggbeeswarm", "0.7.2")'
RUN R -e 'remotes::install_version("egg", "0.4.5")'
RUN R -e 'remotes::install_version("cowplot", "1.1.1")'
RUN R -e 'remotes::install_version("ComplexUpset", "1.3.3")'
RUN R -e 'remotes::install_version("stringr", "1.5.0")'
RUN R -e 'remotes::install_version("reshape2", "1.4.4")'
RUN R -e 'remotes::install_version("ggpubr", "0.6.0")'
RUN R -e 'remotes::install_version("patchwork", "1.2.0")'
RUN R -e 'remotes::install_version("survival", "3.5-5")'
RUN R -e 'remotes::install_version("lubridate", "1.9.2")'
RUN R -e 'remotes::install_version("ggsurvfit", "0.1.0")'

# R Bioconductor packages:
#RUN R -e 'BiocManager::install("BSgenome")'
RUN R -e 'BiocManager::install("GenomicRanges")'
RUN R -e 'BiocManager::install("ggbio", "1.38.0")'
RUN R -e 'BiocManager::install("EnsDb.Hsapiens.v86", "2.99.0")'
