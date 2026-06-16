
#' @export
setMethod("[", "data.struct",
  function(x, i, j, ..., drop = FALSE) { 
    if(missing(i)) i <- seq_len(nrow(x))
    if(missing(j)) j <- seq_len(ncol(x))
    if(is.character(j)) {
      ds_names = colNamesDataStruct(x)
      j <- match(j, ds_names)
    }
    if(drop && length(j) == 1) {
      return( x[[j]][i] )
    }
    i <- as.integer(i) - 1L
    j <- as.integer(j) - 1L
    if(anyNA(i)) stop("Undefined lines")
    if(anyNA(j)) stop("Undefined columns")
    new("data.struct", ptr = extractDataStruct(x, i, j), type = 0L)
  }
)


