#' @export
setGeneric("GRM", function(x, ...) standardGeneric("GRM"))

#' Genetic Relationship Matrix
#' 
#' @description  Compute the Genetic Relationship Matrix (GRM) between given SNPs.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param ... extra arguments for some methods. Could be a file name or 
#' a bool speciifying to use float
#'  
#' @details
#' The GRM is a symmetric square matrix of dimension equal to the number of individuals.
#' Each entry can be interpreted as an estimated kinship coefficient between individuals,
#' although some authors might disagree. Note in particular that some entries will be negative.
#' 
#' @return
#' A matrix or mmatrix of GRM values.
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
