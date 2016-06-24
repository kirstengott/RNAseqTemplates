#!/usr/bin/env Rscript

## configuration

workingdir <- "" ## ?getwd()
control    <- "" ## "control, control, control"
treatment  <- "" ## "treatment, treatment, treatment" 
pval       <- "" ## "0.001"



setwd(workingdir)

## First list the R libraries I want to use
libraries        <- c("DESeq2", "parallel", "dplyr", 'tidyr', 'kirsten')


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




countsfls <- list.files("./counts", full.names = TRUE)


countsls <- lapply(countsfls, function(x){
  if(x == countsfls[1]){
      y           <- data.table::fread(x, sep = "\t", data.table = FALSE)[-1, c(1,2,3)]
      colnames(y) <- c('ID', 'transcript_length', basename(x))
      invisible(y)
  } else {
      y           <- data.table::fread(x, sep = "\t", data.table = FALSE)[-1, 3, drop = FALSE]
      colnames(y) <- basename(x)
      invisible(y)
  }})


counts     <- bind_cols(countsls)
tidyCounts <- gather(counts, Sample, Count, -transcript_length, -ID)




normCounts <- tidyCounts %>% group_by(Sample) %>%
                             mutate(rpkm = (1e9 * as.numeric(Count))/(sum(as.numeric(Count)) * as.numeric(transcript_length))) %>%
                               mutate(rpm = (as.numeric(Count))/(sum(as.numeric(Count))/1e6)) %>%
                                  mutate(rpk = (as.numeric(Count)/(as.numeric(transcript_length)/1000))) %>%
                                      mutate(tpm = rpk/(sum(as.numeric(rpk))/1e6)) %>% ungroup() %>% select(-rpk)










rpkmOut <- normCounts %>% select(ID, transcript_length, Sample, rpkm) %>% spread(Sample, rpkm)
tpmOut  <- normCounts %>% select(ID, transcript_length, Sample, tpm) %>% spread(Sample, tpm)
rpmOut  <- normCounts %>% select(ID, transcript_length, Sample, rpm) %>% spread(Sample, rpm)



WriteTable(x = rpkmOut, file = 'output/rpkm.txt')
WriteTable(x = tpmOut, file = 'output/tpm.txt')
WriteTable(x = rpmOut, file = 'output/rpm.txt')



countsData <- data.frame(row.names = counts$ID, (counts %>% select(-ID, -transcript_length)), check.names = FALSE)

colData           <- data.table::fread("./sample_info", sep = "\t", data.table = FALSE, header = FALSE, col.names = c('condition', 'replicate'))
rownames(colData) <- colData$replicate


countsDataOrdered <- countsData[,rownames(colData)]

## run deseq on anything (for the whole PCA plot)
dds_all <- DESeqDataSetFromMatrix(countData = countsDataOrdered,
                                  colData = colData,
                                  design = ~ condition)





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




dds.res <- lapply(contrast, function(x){
    x1                 <- x[1]
    x2                 <- x[2]
    small.counts.data  <- countsData[,c(grep(x1, colnames(countsData)), grep(x2, colnames(countsData)))]
    small.col.data     <- colData[c(grep(x1, rownames(colData)), grep(x2, rownames(colData))),]
    print(x)
    DESeqDataSetFromMatrix(countData = small.counts.data,
                           colData = small.col.data,
                           design = ~ condition)
})
    


contrast_names <- unlist(lapply(contrast, paste, collapse = "/"))
names(dds.res) <- contrast_names


dds.res <- lapply(contrast_names, function(x){
    print(x)
    dds <- dds.res[[x]]
    DESeq(dds)})

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
})


names(all.results) <- contrast_names


## Make a long and wide results DF

all_results_tidyDF <- bind_rows(all.results)

all_results_DF     <- all_results_tidyDF %>% unite(ID_all, contrastID, dea_ID) %>% spread(key = ID_all, value = dea_Value)


## make PCA info

rld     <- rlogTransformation(dds_all, blind = TRUE)
rld_all <- mclapply(dds.res, rlogTransformation, blind = TRUE)

    
## Save everything we have made to Rdata files

save(file = "workspace_images/counts.Rdata", list = c('counts', 'tidyCounts', 'normCounts'))
save(file = "workspace_images/config.Rdata", list = c('title', 'workingdir', 'control', 'treatment', 'pval'))
save(file = "workspace_images/DESeq2All.Rdata", list = c('title', 'workingdir', 'pval', 'tidyCounts', 'normCounts', 'countsData', 'contrast', 'dds_all', 'all_results_tidyDF', 'all_results_DF', 'contrast_names', 'dds.res', 'rld', 'rld_all'))
