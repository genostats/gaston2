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
GRM <- function(x, ...) UseMethod("GRM")

#' @export
setMethod("GRM", c(x = "snp.matrix"),
  function(x, ..., filename, usefloat = FALSE) {
    if(missing(filename)) {
      K <- grm_(x@ptr) 
    } else {
      if (is.null(usefloat)) {
        K <- grm_mmatrix(x@ptr, filename)
      } else {
        K <- grm_mmatrix(x@ptr, filename, usefloat)
      }
    }
    ids <- ind.stats(x)
    uid <- uniqueid(ids$famid, ids$id)
    rownames(K) <- colnames(K) <- uid
    K
  }
)
