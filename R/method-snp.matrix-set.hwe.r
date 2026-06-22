#' Hardy-Weinberg Equilibrium
#' 
#'
#' @export
set.hwe <- function(x, method, verbose) UseMethod("set.hwe")

#' @export
set.hwe.snp.matrix <- function(x, method = c("chisquare", "yates", "exact"), verbose = FALSE) {
  test <- match.arg.index(method, TRUE)
  set_hwe(x@ptr, test, FALSE)
}

