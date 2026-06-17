#' @export 
setMethod("[<-", "data.struct",
  function(x, i, j, ..., value) {
    dd <- dim(x)
    nr <- dd[1]
    nc <- dd[2]
    if(missing(j)) j <- 1:nc
    # a single colunm is concerned
    if(length(j) == 1) { # x[i, col] <- value 
      if(missing(i)) {
        x[[j]] <- value
      } else {
        # could use 'x[[j]][i] <- value' but this is cleaner + avoid column type change (when cast is possible)
        if(is.character(j)) {
          ds_names <- colNamesDataStruct(x)
          j <- match(j, ds_names)
        }
        j <- as.integer(j) - 1L
        if(is.na(j)) stop("Undefined columns")
        if(j > ncolDataStruct(x)) stop("Column index out of range")
        setCellsDataStruct(x, i - 1L, j, value);
      }
    } else { # several cols
      if(missing(i)) i <- 1:nr
      hasCols <- !is.null(dim(value)) && length(dim(value)) == 2
      for(k in seq_along(j)) {
        if(hasCols) 
          x[[ j[k] ]][i] <- value[, k]
        else
          x[[ j[k] ]][i] <- value
      }
    }
    return(x)
  }
)



