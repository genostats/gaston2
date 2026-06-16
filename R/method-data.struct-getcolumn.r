#' @export
setMethod("$", "data.struct",
  function(x, name) {
    ds_names <- colNamesDataStruct(x)
    if(!(name %in% ds_names))
      NULL
    else
      getColAsSEXPDataStruct(x, name)
  }
)


#' @export
setMethod("[[", "data.struct",
  function(x, i, j, ...) { # j should always be missing
    if (length(i) != 1) stop("Multiple column subsetting not allowed with [[")
    ds_names <- colNamesDataStruct(x)
    if(is.numeric(i) && i <= ncol(x)) { # il faudrait faire une fonction qui prend un numéro de colonne
      return(getColAsSEXPDataStruct(x, ds_names[i]))
    }
    i <- as.character(i)
    if(!(i %in% ds_names))
      NULL
    else 
      getColAsSEXPDataStruct(x, i)
  }
)


