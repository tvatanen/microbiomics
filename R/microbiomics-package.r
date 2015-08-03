#' Microbiomics
#' 
#' Implementation of some functions that may come handy in analysing microbiome data. 
#' Following functions are included:
#' 
#' \code{\link{read_metaphlan_table}} read a MetaPhlAn 2.0 output file as an R data frame
#' \code{\link{read_taxon_table}} read a summarize_taxa.py output file as an R data frame
#' \code{\link{pairwise_spearman}} compute pairwise correlations between columns of two data sets
#' \code{\link{get_significant_values}} get significant components among a matrix of measurements
#' \code{\link{Maaslin.wrapper}} run MaAsLin (\url{https://bitbucket.org/biobakery/maaslin/}) given a 
#' set of taxa and metadata.
#' 
#' @name microbiomics-package
#' @docType package
#' 
#' @import pspearman
#' @import magrittr
#' @import stringr
NULL
