
#' @export
setGeneric("select.snps", function(x, condition) standardGeneric("select.snps"))

#' @export
setMethod("select.snps", "snp.matrix",
  function(x, condition) {
    w <- eval_ds(substitute(condition), snp.stats(x))
    miss <- is.na(w)
    n.miss <- sum(miss)
    if(n.miss > 0) {
      warning(paste(n.miss, "SNP(s) with undefined condition are removed from snp.matrix"))
      w <- w & !miss
    }
    x[, w]
  }
)

