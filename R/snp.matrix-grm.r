#' @export
setGeneric("GRM", function(x) standardGeneric("GRM"))

#' @export
setMethod("GRM", c(x = "snp.matrix"), 
  function(x) {
    return(grm_(x@ptr))
  }
)
