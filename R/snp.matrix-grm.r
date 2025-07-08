#' @export
setGeneric("grm", function(x) standardGeneric("grm"))
#je crois qu'il manquait un setGeneric

#' @export
setMethod("grm", c(x = "snp.matrix"), 
  function(x) {
  return(grm_(x@ptr))
  }
)