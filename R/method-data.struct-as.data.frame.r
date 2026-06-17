#' Converting custom datastruct to R data.frame
#' 
#' @name as.data.frame.data.struct
#'
#' @exportS3Method as.data.frame data.struct
as.data.frame.data.struct <- function(x, ...) DataStructToSEXP(x)

#' @family data.struct
#' @exportMethod coerce
setAs("data.struct", "data.frame", function(from) DataStructToSEXP(x))

