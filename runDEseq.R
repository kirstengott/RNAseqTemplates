#!/usr/bin/env Rscript


## First list the R libraries I want to use
libraries        <- c("DESeq2", "parallel", 'kiRsten', 'tidyverse')



## Now read in the libraries
lapply(libraries, function(x){
  print(x)
  library(x, character.only = TRUE, quietly = TRUE)
})


## Create output directories


new.dirs <- c('figures', 'bin', 'workspace_images', 'output')
lapply(new.dirs, function(x){
  if(!dir.exists(x)){
    dir.create(x)}})




tidyCounts <- data.table::fread('output/countsTidy.txt', data.table = FALSE, sep = "\t")

counts <- spread(tidyCounts, Sample, Count)

WriteTable(x = counts, file = 'output/counts.txt')
WriteTable(x = tidyCounts, file = 'output/countsTidy.txt')

normCounts <- tidyCounts %>% group_by(Sample) %>%
  mutate(rpkm = (1e9 * as.numeric(Count))/(sum(as.numeric(Count)) * as.numeric(transcript_length))) %>%
  mutate(rpm = (as.numeric(Count))/(sum(as.numeric(Count))/1e6)) %>%
  mutate(rpk = (as.numeric(Count)/(as.numeric(transcript_length)/1000))) %>%
  mutate(tpm = rpk/(sum(as.numeric(rpk))/1e6)) %>% 
  ungroup() %>% 
  dplyr::select(-rpk)










rpkmOut <- normCounts %>% dplyr::select(ID, transcript_length, Sample, rpkm) %>% spread(Sample, rpkm)
tpmOut  <- normCounts %>% dplyr::select(ID, transcript_length, Sample, tpm) %>% spread(Sample, tpm)
rpmOut  <- normCounts %>% dplyr::select(ID, transcript_length, Sample, rpm) %>% spread(Sample, rpm)



WriteTable(x = rpkmOut, file = 'output/rpkm.txt')
WriteTable(x = tpmOut, file = 'output/tpm.txt')
WriteTable(x = rpmOut, file = 'output/rpm.txt')





countsData <- data.frame(row.names = counts$ID, 
                         (counts %>% dplyr::select(-ID, -transcript_length)), 
                         check.names = FALSE)

countsData <- countsData[, colnames(countsData)[which(!colnames(countsData) %in% samplesRemove)]]




colData           <- data.table::fread("./sample_info", sep = "\t", data.table = FALSE, header = FALSE, col.names = c('condition', 'replicate'))
rownames(colData) <- colData$replicate

colData <- colData[rownames(colData) %in% colnames(countsData), ]


countsDataOrdered <- countsData[ ,rownames(colData)]

## run deseq on anything (for the whole PCA plot)
dds_all <- DESeqDataSetFromMatrix(countData = countsDataOrdered,
                                  colData = colData,
                                  design = ~ condition)



dds_all_res <- DESeq(dds_all)




## run the rest of the analysis


treatment.init <- gsub(" ", "", unlist(strsplit(treatment, split = ", ")))
control.init   <- gsub(" ", "", unlist(strsplit(control, split = ", ")))

    ## Create contrasts, if there are different numbers of treatment and control, throw an error
if(length(treatment.init) == length(control.init)){
    contrast <- lapply(seq(1, length(treatment.init)), function(x){
        c(treatment.init[x], control.init[x])})
} else {
    message("\nNumber of controls does not equal number of treatments.  Exiting program")
}




dds.res <- mclapply(contrast, function(x){
    x1                 <- x[1]
    x2                 <- x[2]
    small.counts.data  <- countsData[,c(grep(x1, colnames(countsData)), grep(x2, colnames(countsData)))]
    small.col.data     <- colData[c(grep(x1, rownames(colData)), grep(x2, rownames(colData))),]
    print(x)
    DESeqDataSetFromMatrix(countData = small.counts.data,
                           colData = small.col.data,
                           design = ~ condition)
}, mc.cores = cores)
    


contrast_names <- unlist(lapply(contrast, paste, collapse = "/"))
names(dds.res) <- contrast_names


dds.res <- mclapply(contrast_names, function(x){
    print(x)
    dds <- dds.res[[x]]
    DESeq(dds)}, mc.cores = cores)

names(dds.res) <- contrast_names


## Run DESEq and extract results

all.results <- mclapply(contrast_names, function(x){
    print(x)
    dds <- dds.res[[x]]
    l <- strsplit(x, split = "/")[[1]]
    y <- data.frame(results(dds, contrast=c("condition",l[1],l[2]), tidy = TRUE))[,c('row', 'baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj')]
    y_f <- gather(y, key = dea_ID, value = dea_Value, -row) %>% dplyr::rename('Gene' = row)
    y_f$contrastID <- x
    y_f
}, mc.cores = cores)


names(all.results) <- contrast_names



all.results.glm <- mclapply(contrast_names, function(x){
  print(x)
  l <- strsplit(x, split = "/")[[1]]
  y <- data.frame(results(dds_all_res, contrast=c("condition",l[1],l[2]), tidy = TRUE))[,c('row', 'baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj')]
  y_f <- gather(y, key = dea_ID, value = dea_Value, -row) %>% dplyr::rename('Gene' = row)
  y_f$contrastID <- x
  y_f
}, mc.cores = cores)


names(all.results.glm) <- contrast_names



## Make a long and wide results DF

all_results_glmTidy <- bind_rows(all.results.glm)
all_results_tidyDF  <- bind_rows(all.results)
all_results_DF      <- all_results_tidyDF %>% unite(ID_all, contrastID, dea_ID) %>% spread(key = ID_all, value = dea_Value)

## write them out
WriteTable(x = all_results_tidyDF, file = 'output/all_genes_expressionTidy.txt')

lfc_table <- all_results_tidyDF %>%
  mutate(significant = ifelse(dea_ID == 'padj' & dea_Value <= pval, yes = 'yes', no = 'no')) %>% # label all significant genes
  group_by(Gene) %>% mutate(num_sigGroup = length(which(significant == 'yes'))) %>% # label the # of times genes that are significant
  ungroup() %>%
  dplyr::filter(num_sigGroup > 0, dea_ID == 'log2FoldChange') %>% # pull out the fold change for genes significant at least 1 time
  dplyr::select(-significant, -num_sigGroup) %>% unite(colname, contrastID, dea_ID) %>% # make the table wide formatted
  spread(colname, dea_Value)

sign.table <- all_results_tidyDF  %>%
  spread(key = dea_ID, value = dea_Value)  %>%
  filter(Gene %in% lfc_table$Gene) %>%
  group_by(contrastID, Gene) %>%
  mutate(sortby = -log(padj, base = c(10))*sign(log2FoldChange)) %>%
  ungroup() %>% gather(dea_ID, dea_Value, -Gene, -contrastID) %>%
  unite(idAll, contrastID, dea_ID) %>% spread(key = idAll, value = dea_Value) %>%
  dplyr::rename('TranscriptID' = Gene)



#Make a table that contains all of the genes in the results
all.genes.table <- all_results_DF %>% dplyr::rename('TranscriptID' = Gene)


WriteTable(file = paste0("./output/all_genes_expression.txt"), x = all.genes.table)
WriteTable(file = paste0("./output/significant_genes.txt"), x = sign.table)

## make PCA info

rld     <- rlogTransformation(dds_all, blind = TRUE)
rld_all <- mclapply(dds.res, rlogTransformation, blind = TRUE, mc.cores = cores)




## make heatmap inputs


## fix the log Fold Change table so it can be clustered

heatmap.input.table           <- data.frame(lfc_table)
rownames(heatmap.input.table) <- lfc_table$Gene
heatmap.input.table$Gene      <- NULL
heatmap.input.table           <- na.omit(data.matrix(heatmap.input.table))
colnames(heatmap.input.table) <- contrast_names



    
## Save everything we have made to Rdata files

save(file = "workspace_images/counts.Rdata", list = c('counts', 'tidyCounts', 'normCounts'))
save(file = "workspace_images/config.Rdata", list = c('title', 'workingdir', 'control', 'treatment', 'pval'))
save(file = "workspace_images/DESeq2All.Rdata", list = c('title', 'workingdir', 'pval', 'tidyCounts',
                                                         'normCounts', 'countsData', 'contrast', 'dds_all',
                                                         'all_results_tidyDF', 'all_results_DF', 'contrast_names',
                                                         'dds.res', 'rld', 'rld_all', 'all_results_glmTidy', 'heatmap.input.table'))
