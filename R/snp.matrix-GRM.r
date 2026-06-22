#' @export
setGeneric("GRM", function(x, ...) standardGeneric("GRM"))

#' Genetic Relationship Matrix
#' 
#' @description  Compute the Genetic Relationship Matrix (GRM) between given SNPs.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param ... extra arguments for some methods.
#'  
#' @details
#'
#' @return
#' A matrix or mmatrix.
#' 
#' @keywords  Genetic Relationship Matrix
#' 
#' @export
setMethod("GRM", c(x = "snp.matrix"),
  function(x, ...) {
    filename <- list(...)$filename
    usefloat <- list(...)$usefloat
    if(is.null(filename)) 
      grm_(x@ptr) # for now only computing in float when given to houba
    else {
      if (is.null(usefloat)) grm_mmatrix(x@ptr, filename)
      else grm_mmatrix(x@ptr, filename, usefloat)# output houba's mmatrix
    }
  }
)
