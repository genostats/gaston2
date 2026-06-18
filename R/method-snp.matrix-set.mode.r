#' @export
setGeneric("set.mode", function(object, mode) standardGeneric("set.mode"))

#' @export
setMethod("set.mode", "snp.matrix", 
  function(object, mode) {
    mode <- match.arg(mode, c("raw.values", "centered", "standardized.mu.sigma", "standardized.p"))
    mode <- toupper(mode)
    mode <- gsub("\\.", "_", mode)
    if(!isnullptr(object@ptr)) {
      setMode(object@ptr, mode)
    }
    else{
      stop("snp.matrix has a broken ptr !")
    }
  }
)

#' @export
setMethod("set.mode", "dose.matrix", 
  function(object, mode) {
    mode <- match.arg(mode, c("raw.values", "centered", "standardized.mu.sigma", "standardized.p"))
    mode <- toupper(mode)
    mode <- gsub("\\.", "_", mode)
    if(!isnullptr(object@ptr)) {
      setModeDosage(object@ptr, mode)
    }
    else{
      stop("dose.matrix has a broken ptr !")
    }
  }
)