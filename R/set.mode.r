
#' @export
set.mode <- function(x, mode) {
  mode <- match.arg(mode, c("raw.values", "centered", "standardized.mu.sigma", "standardized.p"))
  mode <- toupper(mode)
  mode <- gsub("\\.", "_", mode)
  if (class(x) == "snp.matrix") setMode(x@ptr, mode)
  if (class(x) == "dose.matrix") setModeDosage(x@ptr, mode)
  else stop("set.mode failed to identify the type of x")
}