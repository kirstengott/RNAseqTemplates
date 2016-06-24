## RNAseq Analysis Templates

My goal with these templates is to provide a basic pipeline for reproducible RNAseq analysis.

There are two files:

runDEseq.R
analysis.Rmd

runDEseq.R creates the Rdata input for analysis.Rmd

analysis.Rmd creates several type of plot including:

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
kirsten
knitr
