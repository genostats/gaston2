#' @export
setGeneric("snp.stats", function(x) standardGeneric("snp.stats"))

#' @export
setMethod("snp.stats", "snp.matrix", function(x) {
  exportSNPStats(x@ptr)
  new("data.struct", ptr = x@ptr, type = 2L)
})


#' @export
setGeneric("snp.stats<-", function(x, value) standardGeneric("snp.stats<-"))

#' @export
setReplaceMethod("snp.stats", c("snp.matrix", "data.struct"), function(x, value) {
  if(!identical(x@ptr, value@ptr)) 
    stop("Replacing snp.stats(x) by another data.struct is not allowed")
  x
})

