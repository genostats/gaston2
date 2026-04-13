#' Converting custom datastruct to R data.frame
#' 
#' @name as.data.frame.data.struct
#'
#' @exportS3Method as.data.frame data.struct
as.data.frame.data.struct <- function(x, ...) {
  DataStructToDataFrame_(x@ptr)
}

#' @family data.struct
setAs("data.struct", "data.frame", function(from) as.data.frame.data.struct(from))

