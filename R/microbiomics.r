#' A function to read MetaPhlAn 2.0 output file
#' 
#' This functions reads a MetaPhlAn 2.0 output file as an R data frame
#' 
#' @param filename MetaPhlAn output file to be read
#' @param kingdom to read: One of "k__Bacteria", "k__Viruses", "k__Eukaryota"
#' @param lvl what level of taxonomy to read (default = 7, species)
#' 
#' @return a data frame with the MetaPhlAn results
#' 
#' @author Tommi Vatanen <tommivat@@gmail.com>
#' @export
read_metaphlan_table <- function(filename, kingdom = "k__Bacteria", lvl = 7, normalize= TRUE) {
  lvl_identifiers <- list("k__","p__","c__","o__","f __","g__","s__")
  if (!(kingdom %in% c("k__Bacteria", "k__Viruses", "k__Eukaryota"))) {
    stop("Kingdom should be on of the following: k__Bacteria [default], k__Viruses, k__Eukaryota")
  }
  if (!(lvl %in% 2:8)) {
    stop("lvl should be integer between 2 (phylotype) and 8 (subspecies / markers)")
  }
  d <- read.table(filename, header=T)
  rownames(d) <- as.character(d[,1])
  d <- d[,-1]
  d <- d[grep(kingdom, rownames(d)),]
  rows <- which(sapply(rownames(d), function(x) length(strsplit(x,'|',fixed=T)[[1]]))==lvl)
  d <- d[rows,]
  d <- data.frame(t(d))
  if (normalize) {
    d <- d / apply(d,1,sum)
  }
  return(d)
}


#' A function to read taxon table produced by summarize_taxa.py
#' 
#' This functions reads a summarize_taxa.py output file as an R data frame
#' 
#' @param filename taxon table file to read
#' @param lvl taxonomic level to read (default = 6, genera)
#' 
#' @return a data frame with the taxonomic data
#' 
#' @author Tommi Vatanen <tommivat@@gmail.com>
#' @export
read_taxon_table <- function(filename, lvl=6) {
  data <- read.table(filename, header=T, sep="\t", stringsAsFactors = F)
  # add rownames (taxonomy + OTU ID)
  
  rownames(data) <- data[,1]
  data <- data[,-1]
  colnames(data) <- substr(colnames(data),2,100)
  data <- as.data.frame(t(data))
  return(data)
}

#' A function to compute pairwise correlations between two data sets
#' 
#' This function computes pairwise correlations between columns of two data sets.
#' It assumes that the data sets have same samples (on rows) in same order.
#' 
#' @param x data set 1 (samples x features1)
#' @param y data set 2 (samples x features2)
#' 
#' @return a list with following elements
#' \itemize{
#'     \item  the clustering of rows as \code{\link{hclust}} object 
#'     \item r features1 x features2 matrix of pairwise correlations
#'     \item pval features1 x features2 matrix of nominal p-values
#' }
#' 
#' @author Tommi Vatanen <tommivat@@gmail.com>
#' @export
pairwise_spearman <- function(x, y) {
  pvalues <- array(0,c(dim(x)[2],dim(y)[2]))
  spearman <- array(0,c(dim(x)[2],dim(y)[2]))
  
  for (i in 1:dim(x)[2]) {
    for (j in 1:dim(y)[2]) {
      tmp <- spearman.test(x[,i],y[,j], alternative = "t", approximation = "AS89")
      pvalues[i,j] <- tmp$p.value
      spearman[i,j] <- tmp$estimate
    }
  } 
  return(list(r=spearman,pval=pvalues))
}

