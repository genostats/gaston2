#' Converting snp matrices to R matrices
#' 
#' @name as.matrix.snp.matrix
#'
#' @exportS3Method as.matrix snp.matrix
as.matrix.snp.matrix <- function(x, ...) {
  mode <- get.mode(x)
  if(mode == "raw.values")
    return(SNPMatrixToIntegerMatrix(x@ptr))
  else
    return(SNPMatrixToNumericMatrix(x@ptr))
}

#' @family snp.matrix
setAs("snp.matrix", "matrix", function(from) as.matrix.snp.matrix(from))

