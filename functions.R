splitstr2vec <- function(x, split = ' '){
  strsplit(x, split = split)[[1]]
}

prettyLinkFile <- function(x, y, sep = '\n\n'){
  cat(paste0('[', x, ']', '(',  y, ')'), sep)
}

se <- function(x) sqrt(var(x)/length(x))

prettyFilePath <- function(x, prepend = '', sep = '\n\n'){
  cat(paste0(prepend, x), sep = sep)
}

make_pheatmap_groups <- function(hclust_object, numClusters = 4){
  ## Takes and object from 'hclust' and the number of Clusters to assign groups to.
  ## Returns a data frame of the properly formated annotation to use with pheatmap.
  mycl                 <- cutree(hclust_object, k = numClusters)
  cluster.letters      <- LETTERS[seq( from = 1, to = numClusters)]
  clusters             <- paste0("Cluster ", cluster.letters)
  mycols               <- clusters[as.vector(mycl)]
  names(mycols)        <- names(mycl)
  annotation           <- data.frame(Cluster=mycols)
  rownames(annotation) <- names(mycl)
  annotation$Cluster   <- as.factor(annotation$Cluster)
  invisible(annotation)
}

make_pheatmap_df <- function(x, min = -3, max = 3){
  ## Makes a dataframe that takes kindly to pheatmap's color mapping restrctions.
  ## x: a dataframe to make a heatmap from.
  x[x > max] <- max
  x[x < min] <- min
  invisible(x)
}

