#' @export
setGeneric("reset.stats", function(x) standardGeneric("reset.stats"))

#' @export
setMethod("reset.stats", "snp.matrix", function(x) {
  # call with force = TRUE
  computeIndStats(x@ptr, TRUE)
  exportSNPStats(x@ptr, TRUE)
  # recompute HWE if needed
  nn <- names(snp.stats(x))
  if("hwe.chi2"  %in% nn) set.hwe(x, "chi2")
  if("hwe.yates" %in% nn) set.hwe(x, "yates")
  if("hwe.exact" %in% nn) set.hwe(x, "exact")
  invisible(NULL)
})

