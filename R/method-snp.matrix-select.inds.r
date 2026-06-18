
#' @export
setGeneric("select.inds", function(x, condition) standardGeneric("select.inds"))

#' @export
setMethod("select.inds", "snp.matrix",
  function(x, condition) {
    cdt <- substitute(condition)
    # get a data frame with only the relevant columns
    df <- as.data.frame(ind.stats(x), cols = all.names(cdt))
    w <- eval(cdt, df, parent.frame())
    miss <- is.na(w)
    n.miss <- sum(miss)
    if(n.miss > 0) {
      warning(paste(n.miss, "individual(s) with undefined condition are removed from snp.matrix"))
      w <- w & !miss
    }
    x[w, ]
  }
)

