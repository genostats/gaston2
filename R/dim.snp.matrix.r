#' @exportS3Method dim snp.matrix
dim.snp.matrix <- function(x) dimSNPmatrix(x@ptr)
