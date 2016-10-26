#' Microbiomics
#' 
#' Implementation of some functions that may come handy in analysing microbiome data. 
#' 
#' Following functions are included:
#' \describe{
#'  \item{\code{\link{read_metaphlan_table}}}{read a MetaPhlAn 2.0 output file as an R data frame}
#'  \item{\code{\link{jaccard}}}{compute Jaccard index between given microbial profiles}
#'  \item{\code{\link{read_taxon_table}}}{read a summarize_taxa.py output file as an R data frame}
#'  \item{\code{\link{pairwise_spearman}}}{compute pairwise correlations between columns of two data sets}
#'  \item{\code{\link{get_significant_values}}}{get significant components among a matrix of measurements}
#'  \item{\code{\link{Maaslin.wrapper}}}{run MaAsLin (\url{https://bitbucket.org/biobakery/maaslin/}) given a
#' set of taxa and metadata}
#' }
#' 
#' @import pspearman
#' @import magrittr
#' @import stringr
#' 
#' @docType package
#' @name microbiomics-package
NULL
#> NULL