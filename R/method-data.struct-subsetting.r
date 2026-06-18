
#' @export
setMethod("[", "data.struct",
  function(x, i, j, ..., drop = TRUE) { browser();
    if(missing(i)) i <- seq_len(nrow(x))
    if(missing(j)) j <- seq_len(ncol(x))
    if(is.character(j)) {
      ds_names <- colNamesDataStruct(x)
      j <- match(j, ds_names)
    }
    if(any(i <= 0L)) i <- seq_len(nrow(x))[i]
    if(any(j <= 0L)) j <- seq_len(ncol(x))[j]
    i <- as.integer(i) - 1L
    j <- as.integer(j) - 1L
    if(anyNA(i)) stop("Undefined lines")
    if(anyNA(j)) stop("Undefined columns")
    r <- new("data.struct", ptr = extractDataStruct(x, i, j), type = 0L)
    if(drop && length(j) == 1) 
      r[[1]] 
    else
      r
  }
)


