
#' @exportMethod coerce
setAs("data.frame", "data.struct", function(from) new("data.struct", ptr = DataFrameToDataStructR_(from), type = 0L))

# S3 or S4?
#' @export
as.data.struct <- function(x) UseMethod("as.data.struct")

#' @export
as.data.struct.data.frame <- function(x) new("data.struct", ptr = DataFrameToDataStructR_(x), type = 0L)



