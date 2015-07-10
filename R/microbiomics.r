#' A function to read MetaPhlAn 2.0 output file
#' 
#' This functions read a MetaPhlAn 2.0 output file as an R data frame
#' 
#' @param filename MetaPhlAn output file to be read
#' @param kingdom to read: One of "k__Bacteria", "k__Viruses", "k__Eukaryota"
#' @param lvl what level of taxonomy to read (default = 7, species)
#' 
#' @return a data frame with the MetaPhlAn results
#' 
#' @author Tommi Vatanen tommivat@gmail.com
#' @export
read_metaphlan_table <- function(filename, kingdom="k__Bacteria", lvl=7) {
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
  return(data.frame(t(d)))
}