
#' @export
setGeneric("select.inds", function(x, condition) standardGeneric("select.inds"))

#' @export
setMethod("select.inds", "snp.matrix",
  function(x, condition) {
    w <- eval_ds(substitute(condition), ind.stats(x))
    miss <- is.na(w)
    n.miss <- sum(miss)
    if(n.miss > 0) {
      warning(paste(n.miss, "individual(s) with undefined condition are removed from snp.matrix"))
      w <- w & !miss
    }
    x[w, ]
  }
)

