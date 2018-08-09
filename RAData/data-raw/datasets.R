library(knitr)
knit2html('collateTables.Rmd',envir=topenv(),output='../inst/extdata/Logfiles')
ls()
keepDataObjects(c('radup', 'radedup'))

#' RAData
#' A data package for study RAData
#' @docType package
#' @aliases RAData-package
#' @title RAData
#' @name RAData
#' @description a description of the package.
#' @details Additional details.
#' @import data.table
#' @import stringr
#' @import knitr
#' @seealso \link{mydataset}
NULL

#' Data from single cell RNAseq of RA patients, before UMI deduplications
#'@name radup
#'@docType data
#'@title Bulk and single cell data, before deduplication
#'@format a list of length three containing
#'\describe{
#'\item{exprs}{matrix of log2 CPM expression, rows are cells, columns are genes}
#'\item{cData}{data.table of cell data}
#' \item{fData}{data.table of feature data}
#'}
#'@source Dan Lui @ Robinson Lab, Stanford.
#'@seealso \link{RAData}
#' @aliases radedup

