
#' @export
get.mode <- function(x) {
  if (class(x) == "snp.matrix") {
    mode <- getMode(x@ptr)
  }
  else {
    mode <- getModeDosage(x@ptr)
  }
  mode <- tolower(mode)
  gsub("_", ".", mode)
}
