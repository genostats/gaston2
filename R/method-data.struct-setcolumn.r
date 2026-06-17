#' @export
setMethod("$<-", "data.struct",
 function(x, name, value) {
   value <- as.vector(value)
   # on mime le comportement d'un data frame
   nr <- nrowDataStruct(x)
   if(length(value) != nr) {
     if(length(value) == 1 && nr > 0) {
       value <- rep(value, nr)
     } else { # on autorise à ajouter une colonne de taille quelconque quand la dim est 0x0 !
       if(nr > 0 || ncolDataStruct(x) > 0) stop("Length of replacement is ", length(value), ", but data.struct has ", nr, " rows")
     }
   }
   addColDataStruct(x, name, value)
   return(x)
  }
)

#' @export 
setMethod("[[<-", "data.struct",
  function(x, i, j, ..., value) {
    if(missing(j)) { # i is column index, call is x[[i]] <- value (col replacement)
      # check if i is ok
      if(length(i) != 1) stop("Can only set one column at a time")
      if(is.numeric(i)) { # i est numeric
        if(i > ncol(x) + 1) stop("Can't assign to column ", i, " in a data.struct with ", ncol(x), " columns")
        # check length of value
        if(length(value) != nrow(x)) {
          if(length(value) == 1) 
            value <- rep(value, nrow(x))
          else
            stop("Length of replacement is ", length(value), ", but data.struct has ", nrow(x), " rows")
        }
        # get or construct name
        if(i > ncol(x)) 
          nn <- sprintf("V%d", i)
        else
          nn <- colNamesDataStruct(x)[i]
      } else {
        nn <- as.character(i)
      }
      # go for it
      addColDataStruct(x, nn, value)
    } else {
      # call v[[i, j]] <- value  (cell replacement)
      if(length(i) != 1 || length(j) != 1) stop("Only one value should be replaced")
      if(is.character(j)) {
        ds_names <- colNamesDataStruct(x)
        j <- match(j, ds_names)
      }
      if(is.na(i) || is.na(j)) stop("Undefined cell coordinates")
      if(i <= 0L || j <= 0L) stop("Coordinates should be positive")
      setCellDataStruct(x, i - 1L, j - 1L, value)
    }
    x
  }
)



