#' @export
setGeneric("snp.stats", function(x) standardGeneric("snp.stats"))

#' @export
setMethod("snp.stats", "snp.matrix", function(x) {
  new("data.struct", ptr = x@ptr, type = 2L)
})

