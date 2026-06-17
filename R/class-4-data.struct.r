#' @name data.struct
#' @rdname data.struct
#' @title Data Structure
#' 
#' @description A 'data.frame' like class of object for SNP stats etc.
#' 
#' @details 
#' @exportClass data.struct
setClass("data.struct", slots = c(ptr = "externalptr", type = "integer")) # equivalent to nullptr
# type = 0 pour des data struct autonomes,
#        1 pour les stats des individus d'une snp matrix
#        2 pour les stats des SNPs d'une snp matrix

setMethod("show", "data.struct", 
  function(object) {
    if(isnullptr(object@ptr)) {
      cat("A data.struct with a broken external ptr\n")
    } else {
      dd <- dimDataStruct(object)
      cat("A", dd[1], "x", dd[2], "data.struct\n")
      # list all Columns present
      vecNames <- colNamesDataStruct(object)
      for (name in vecNames) {  # for every col
        cat("$ ", name, ": ")
        type <- colTypeDataStruct(object, name)
        cat(tolower(type), " ")
        chartaken = 5 + nchar(name) + nchar(type)
        # also cat beginning of columns... with a threshold of end of line 
        # je voulais aller que aux 3/4 mais peut etre que c'est nul en fait
        maxcharleft <- signif((getOption("width") * 0.75) - chartaken, digits=2) - 1 # to have room for \n 
          cat(showColDataStruct(object, name, maxcharleft))
        cat("\n")
      }
    }
    # TODO : maybe add a way to show the parent matrix if it exists ? 
  } 
)


#' @name data.struct
#' @description A constructor for data.struct from a dataframe
#' @export
data.struct <- function(...) {
  df <- data.frame(...)
  new("data.struct", ptr = DataFrameToDataStructR_(df), type = 0L)
}
