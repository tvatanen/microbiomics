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
#'     \item r: features1 x features2 matrix of pairwise correlations
#'     \item p: features1 x features2 matrix of nominal p-values
#'     \item q: features1 x features2 matrix of fdr corrected p-values (q-values)
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
  qvalues <- array(p.adjust(pvalues), dim(spearman))
  return(list(r=spearman,p=pvalues,q=qvalues))
}

#' A function to select significant components among matrix of measurements
#' 
#' This function takes a matrix of values (e.g. correlations from pairwise spearman)
#' and corresponding p-values and/or q-values. Given the threshold for p- and q-values 
#' it will return reduced matrices of significant values where each row and column 
#' contain at least one significant value (while retaining the original row and 
#' column names).
#' 
#' @param values matrix of values of interest (e.g. correlations)
#' @param p_values matrix of p-values for the matrix in the first argument
#' @param q_values (optional) matrix of q-values for the matrix in the first argument. If q-values are given, no p-value thresholding is done.
#' @pvalue_threshold threshold for significance for the p-values (default 0.01)
#' @qvalue_threshold threshold for significance for the q-value (default 0.1)
#' 
#' @return a list with following elements
#' \itemize{
#'     \item values: matrix of significant values in the matrix of the first argument 
#'     \item p: matrix of p-values for the values above
#'     \item q: matrix of q-values for the values above (if q-values are given as argument)
#' }
#' 
#' @author Tommi Vatanen <tommivat@@gmail.com>
#' @export
get_significant_values <- function(values, p_values, q_values=NA, pvalue_threshold=0.01, qvalue_threshold=0.1) {
  
  if (all(is.na(q_values))) {
    ind <- which(p_values < pvalue_threshold, arr.ind = T)
    
    values_new = values[unique(ind[,1]),unique(ind[,2])]
    pvalues_new = p_values[unique(ind[,1]),unique(ind[,2])]
    return(list(values=values_new,p_values=pvalues_new))
    
  } else {
    ind <- which(q_values < qvalue_threshold, arr.ind = T)
    
    values_new = values[unique(ind[,1]),unique(ind[,2])]
    pvalues_new = p_values[unique(ind[,1]),unique(ind[,2])]
    qvalues_new = q_values[unique(ind[,1]),unique(ind[,2])]
    
    return(list(values=values_new,q_values=qvalues_new,p_values=pvalues_new))
    
  }
}
