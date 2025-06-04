
#' @export
get.mode <- function(x) {
  mode <- getMode(x@ptr)
  mode <- tolower(mode)
  gsub("_", ".", mode)
}
