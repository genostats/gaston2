#' LD thinning
#' 
#'
#' @export
LD.thin <- function(x, threshold, max.dist, beg, end, which.snps, dist.unit, extract, keep, ...) UseMethod("LD.thin")

#' @export
LD.thin.snp.matrix <- function(x, threshold, max.dist, beg = 1L, end = ncol(x), which.snps, dist.unit = c("bases", "cM"), extract = TRUE, 
                               keep = c("right", "left", "random", "priority"), ..., priority) {
  if(missing(which.snps)) {
    which.snps <- rep(FALSE, ncol(x))
    which.snps[beg:end] <- TRUE
  }
  which.snps <- as.logical(which.snps)
  if(length(which.snps) != ncol(x)) 
    stop("'which.snps' should be of length ncol(x)")

  # remove SNPs with no data and monomorphic SNPs
  ss <- snp.stats(x)
  which.snps <- which.snps & (ss$callrate > 0) & (ss$maf > 0)

  dist.unit <- match.arg(dist.unit)
  if(dist.unit == "bases") {
    if(missing(max.dist)) max_dist_bp <- 5e5 else max_dist_bp <- max.dist
    max_dist_cm <- 0
  } else {
    check.dist(x)
    max_dist_bp <- 0
    if(missing(max.dist)) max_dist_cm <- 0.5 else max_dist_cm <- max.dist
  }
  if(max_dist_cm > 100) 
    warning("max.dist value seems very high for dist.unit = \"cM\"")

  keep <- match.arg(keep)
  if(keep == "left") {
    LD_thin_left(x@ptr, threshold, max_dist_bp, max_dist_cm, which.snps)
  } else if(keep == "right") {
    LD_thin_right(x@ptr, threshold, max_dist_bp, max_dist_cm, which.snps)
  } else if(keep == "random") {
    LD_thin_random(x@ptr, threshold, max_dist_bp, max_dist_cm, which.snps)
  } else {
    priority <- as.numeric(priority)
    if(length(priority) != ncol(x))
      stop("'priority' should have length equal to ncol(x)")
    LD_thin_priority(x@ptr, threshold, max_dist_bp, max_dist_cm, which.snps, priority)
  }

  if(!extract) 
    which.snps
  else
    x[,which.snps]
}

