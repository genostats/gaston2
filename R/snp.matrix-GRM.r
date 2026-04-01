#' @export
setGeneric("GRM", function(x, filename = NULL) standardGeneric("GRM"))

#' Genetic Relationship Matrix
#' 
#' @description  Compute the Genetic Relationship Matrix (GRM) between given SNPs.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param filename if non-missing, the GRM matrix will be a mmatrix object
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
setMethod("GRM", c(x = "snp.matrix", filename = NULL),
  function(x, filename) {# no files, output R matrix
    return(grm_(x@ptr))
  }
)

#' @export
setMethod("GRM", c(x = "snp.matrix", filename = "character"),
  function(x, filename) {# output houba's mmatrix
    return(grm_mmatrix(x@ptr, filename))
  }
)
