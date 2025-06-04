
#' @export
set.mode <- function(x, mode) {
  mode <- match.arg(mode, c("raw.values", "centered", "standardized.mu.sigma", "standardized.p"))
  mode <- toupper(mode)
  mode <- gsub("\\.", "_", mode)
  setMode(x@ptr, mode)
}
