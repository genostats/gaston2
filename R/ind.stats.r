#' @export
setGeneric("ind.stats", function(x) standardGeneric("ind.stats"))

#' @export
setMethod("ind.stats", "snp.matrix", function(x) {
  computeIndStats(x@ptr)
  new("data.struct", ptr = x@ptr, type = 1L)
})

