Codes to generate methylation plots
================
Vahid Akbari

This R markdown includes codes to generate methylation related figures,
inlcuding allelic and non-allelic methylation.

# Import required libraries

``` r
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(M3C))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(patchwork))
```

# Import data and plot overal methylation

This will generate overall methylation in tumour samples and normal WGBS
cases.

![](Methylation_Figures_GitHub_files/figure-gfm/Overall_Methylation-1.png)<!-- -->

Overall methylation with respect to mutation status in TET and IDH at
CGIs
![](Methylation_Figures_GitHub_files/figure-gfm/TET_IDH_Mutation_Methylation_CGIs-1.png)<!-- -->

# t-SNE based on tumour type

![](Methylation_Figures_GitHub_files/figure-gfm/tSNE_Tumour_Types-1.png)<!-- -->

# t-SNE based on biopsy site

![](Methylation_Figures_GitHub_files/figure-gfm/tSNE_Biopsy_Sites-1.png)<!-- -->

# Fraction of allelic DMRs based on copy number status

![](Methylation_Figures_GitHub_files/figure-gfm/aDMRs_Count_Genomic_Regions-1.png)<!-- -->

# RET methylation plot

    ## Rows: 1515522 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (8): ID, Sample, Lib, Type, Class, Allele, Chromosome, Gene
    ## dbl (2): Position, Methylation
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'ID', 'Type', 'Class', 'Allele', 'Position'. You can override using the `.groups` argument.
    ## `geom_smooth()` using formula = 'y ~ x'

![](Methylation_Figures_GitHub_files/figure-gfm/RET_Methylation-1.png)<!-- -->

# CDKN2A methylation plot

    ## `geom_smooth()` using formula = 'y ~ x'

![](Methylation_Figures_GitHub_files/figure-gfm/CDKN2A_Methylation-1.png)<!-- -->

# RET Expression box plots

    ## Rows: 35138 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): ID, Sample, Gene, Type, Category
    ## dbl (1): TPM
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'ID', 'Gene', 'Type'. You can override using the `.groups` argument.

![](Methylation_Figures_GitHub_files/figure-gfm/RET_Expression-1.png)<!-- -->

# CDKN2A Expression box plots

![](Methylation_Figures_GitHub_files/figure-gfm/CDKN2A_Expression-1.png)<!-- -->
