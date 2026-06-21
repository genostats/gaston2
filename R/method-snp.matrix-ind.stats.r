#' @export
setGeneric("ind.stats", function(x) standardGeneric("ind.stats"))

#' @export
setMethod("ind.stats", "snp.matrix", function(x) {
  computeIndStats(x@ptr)
  new("data.struct", ptr = x@ptr, type = 1L)
})


#' @export
setGeneric("ind.stats<-", function(x, value) standardGeneric("ind.stats<-"))

#' @export
setReplaceMethod("ind.stats", c("snp.matrix", "data.struct"), function(x, value) {
  if(!identical(x@ptr, value@ptr)) 
    stop("Replacing ind.stats(x) by another data.struct is not allowed")
  x
})

# all the work is done by R when it constructs 'value'. 
# A call like ind.stats(x)$col <- b is transformed by
# ind.stats<-(x, [the data struct obtained by tmp <- ind.stats(x); z$col <- b])

# this allows to add/modify columns in ind.stats but not to replace
# it with a whole other data struct. I think to it as a feature, we could change this.
