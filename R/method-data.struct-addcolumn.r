#' @export
setGeneric("addcolumn", function(x, ...) standardGeneric("addcolumn"))

#' Converting a vector with a name into a column in the datastruct
#' @name addcolumn
#'
setMethod("addcolumn", "data.struct",
  function(x, ...) {
    args <- list(...)
    df_idx <- which(sapply(args, is.data.frame))
    if (length(df_idx) > 1) stop("Please provide only one data.frame")
    if (length(df_idx) == 1) { # I want only one dataframe
      dataframe <- args[[df_idx[1]]] 
      # TODO : check if a name/list of names was given
      i <- 1
      # j'itère sur les colonnes
      while (i <= length(dataframe)) {
      addcolDataStruct(x@ptr, names(dataframe)[i], dataframe[[i]])
      i <- i + 1
      }
    } else {
    colname_idx <- which(sapply(args, is.character))
    vec <- args$vec # is.vector not working well

    if (length(colname_idx) == 0) stop("Please name your new column")
    if (is.null(vec)) stop("Please give values to link to your column (using vec=...)")
    addcolDataStruct(x@ptr, args[[colname_idx[1]]], vec)
    }
  }
)