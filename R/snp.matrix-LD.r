#' Linkage Disequilibrium
#' 
#' @description  Compute Linkage Disequilibrium (LD) between given SNPs.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param lim  Range of SNPs for which the LD is computed
#' @param lim2  (Optional) Second range of SNPs (see Details)
#' @param measure  The LD measure
#' @param filename if non-missing, the LD matrix will be a mmatrix object 
#' @param float if TRUE, computations are done in float (only if filename is provided)
#'  
#' @details
#' If \code{lim2} is missing, the LD is computed between all SNPs with indices between \code{lim[1]} and \code{lim[2]};
#' else, the LD is computed between the SNPs in the range given by \code{lim} and those in the range given by \code{lim2}.
#' 
#' The default is to use LD moment estimators (which are faster but less precise than maximum likelihood estimators
#' obtained with the EM algorithm).
#' 
#' @return
#' A matrix of LD values.
#' @seealso  \code{\link{LD.thin}},  \code{\link{LD.plot}}
#' 
#' @keywords  Linkage Disequilibrium
#' 
#' @export 
setGeneric("LD", function(x, lim, lim2, measure, method, ...) standardGeneric("LD"))

#' @export 
setMethod("LD", c(x = "snp.matrix"), 
  function(x, lim, lim2, measure = c("r2", "r", "D"), method = c("moments", "EM"), filename, float) {

    lim <- as.integer(lim)
    if(length(lim) == 1) lim = c(lim, lim)
    if(length(lim) != 2) stop('lim must be a vector of length 2')
    if(any(lim < 1) | any(lim > ncol(x)) | lim[1] > lim[2]) stop("inappropriate range in lim")
    lim <- lim - 1L

    if(!missing(lim2)) {
      lim2 <- as.integer(lim2)
      if(length(lim2) == 1) lim2 = c(lim2, lim2)
      if(length(lim2) != 2) stop('lim2 must be a vector of length 2')
      if(any(lim2 < 1) | any(lim2 > ncol(x)) | lim2[1] > lim2[2]) stop("inappropriate range in lim2")
      lim2 <- lim2 - 1L
    }

    method <- match.arg(method)
    measure <- match.arg(measure) 
   
    if(measure == "r" || measure == "r2")  
      r.scale <- TRUE
    else
      r.scale <- FALSE
 
    if(missing(lim2) || all(lim2 == lim)) {  
      if(missing(filename)) {
        if(method == "moments") 
          A <- LD_square_moments(x@ptr, lim[1], lim[2], r.scale)
        else
          A <- LD_square_EM(x@ptr, lim[1], lim[2], r.scale)
      } else {
        if(method == "moments") 
          A <- LD_square_moments_mmatrix(x@ptr, lim[1], lim[2], r.scale, float)
        else
          A <- LD_square_EM_mmatrix(x@ptr, lim[1], lim[2], r.scale, float)
      }
    } else {
      if(missing(filename)) {
        if(method == "moments") 
          A <- LD_chunk_moments(x@ptr, lim[1], lim[2], lim2[1], lim2[2], r.scale)
        else
          A <- LD_chunk_EM(x@ptr, lim[1], lim[2], lim2[1], lim2[2], r.scale)
      } else {
        if(method == "moments") 
          A <- LD_chunk_moments_mmatrix(x@ptr, lim[1], lim[2], lim2[1], lim2[2], r.scale, float)
        else
          A <- LD_chunk_EM_mmatrix(x@ptr, lim[1], lim[2], lim2[1], lim2[1], r.scale, float)
      }
    }
  
    if(measure == "r" || measure == "D") return(A)

    if(measure == "r2") {
      if(is(A, "mmatrix")) {
        houba::inplace.prod(A, A)
      } else {
        A <- A**2
      }
      return(A)
    }
    stop("Wut?")
  }
)

