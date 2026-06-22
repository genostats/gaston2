#' @export
setGeneric("reset.stats", function(x) standardGeneric("reset.stats"))

#' @export
setMethod("reset.stats", "snp.matrix", function(x) {
  # call with force = TRUE
  computeIndStats(x@ptr, TRUE)
  exportSNPStats(x@ptr, TRUE)
  # recompute HWE if needed
  if("hwe.chi2" %in% names(snp.stats(x))) set.hwe(x, "chisquare")
  if("hwe.yates" %in% names(snp.stats(x))) set.hwe(x, "yates")
  if("hwe.exact" %in% names(snp.stats(x))) set.hwe(x, "exact")
  invisible(NULL)
})

