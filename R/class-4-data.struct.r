#' @name data.struct
#' @rdname data.struct
#' @title Data Structure
#' 
#' @description A container for columns with dtaas
#' 
#' @details the matrixptr member is a way to stop R's GC from deleting the SNPMatrix that possibly owns the DataStruct
#' (for ex in snp.stats(...TRUE) or ind.stats(...TRUE))
#' @exportClass data.struct
setClass("data.struct", slots = c(ptr = "externalptr", matrixptr = "externalptr"), prototype=c(matrixptr=new("externalptr"))) #equivalent to nullptr

setMethod("show", "data.struct", 
  function(object) {
    if(isnullptr(object@ptr)) {
      cat("A data.struct with a broken external ptr\n")
    } else {
      cat("A data.struct with", ncolDataStruct(object@ptr), "columns\n")
      # list all Columns present
      vecNames <- colNamesDataStruct(object@ptr)
      for (name in vecNames) {  # for every col
        cat("$ ", name, ": ")
        type <- getcolTypeDataStruct(object@ptr, name)
        cat(tolower(type), " ")
        chartaken = 5 + nchar(name) + nchar(type)
        # also cat beginning of columns... with a threshold of end of line 
        # je voulais aller que aux 3/4 mais peut etre que c'est nul en fait
        maxcharleft <- signif((getOption("width") * 0.75) - chartaken, digits=2) - 1 # to have room for \n 
          cat(showcolDataStruct(object@ptr, name, maxcharleft))
        cat("\n")
      }
    }
    # TODO : maybe add a way to show the parent matrix if it exists ? 
  } 
)


#' @name data.struct
#' @description A constructor for data.struct from a dataframe
#' @export
data.struct <- function(df) { # could be "..." later on to feed it something else
  if (!is.data.frame(df)) {
    stop("Please feed it a valid data.frame")
  } else  {
    new("data.struct", ptr = DataFrameToDataStructR_(df))
  }
}

# )
# TO THINK : est-ce que je veux un lien avec houba ? Parce qu'a priori non ms 
# peut etre si les stats sont très grosses ? Utiliser le threshold de houba ?