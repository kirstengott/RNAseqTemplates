## RNAseq Analysis Templates

My goal with these templates is to provide a basic pipeline template for reproducible analysis.

Files:

vanillaAnalysisTemplate.Rmd - basic .Rmd outline to help keep me accountable

runDEseq.R - creates the Rdata input for analysis.Rmd

RNAseqAnalysis.Rmd - creates several type of plot including:

-PCA plot

-MA plot

-Volcano plot

-Heatmaps


Libraries Used:

ggplot2

DESeq2

matrixStats

gplots

RColorBrewer

parallel

lattice

plyr

gtools

dplyr

tidyr

kiRsten (available [here](https://github.com/kirstengott/kiRsten) )

knitr
