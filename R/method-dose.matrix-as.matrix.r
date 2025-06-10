#' Converting dose matrices to R matrices
#'
#' @name as.matrix.dose.matrix
#'
#' @exportS3Method as.matrix dose.matrix
as.matrix.dose.matrix <- function(x, ...) {
  mode <- get.mode(x)
  # if(mode == "raw.values") return(DoseMatrixToNumericMatrix(x@ptr))
  # TODO : see when is it gonna send back integers to me ? 
  # maybe when harcalled ? what would be the mode for that ? 
  # for now will only return in Numeric
  return(DoseMatrixToNumericMatrix(x@ptr))
}

#' @family dose.matrix
setAs("dose.matrix", "matrix", function(from) as.matrix.dose.matrix(from))
