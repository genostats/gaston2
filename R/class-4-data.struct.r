#' @name data.struct
#' @rdname data.struct
#' @title Data Structure
#' 
#' @description A 'data.frame' like class of object for SNP stats etc.
#' 
#' @details Methods have been defined for usual operations.
#'
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
      short_type <- c(INT = "int", DOUBLE = "dbl", FLOAT = "ftl", STRING = "str", BOOL = "boo")
      dd <- dimDataStruct(object)
      cat("A", dd[1], "x", dd[2], "data.struct\n")
      # list all Columns present
      vecNames <- colNamesDataStruct(object)
      le.vn <- max(nchar(vecNames))
      for (name in vecNames) {  # for every col
        cat(sprintf(" $ %*s: ", le.vn, name))
        type <- colTypeDataStruct(object, name)
        cat(short_type[type])
        chartaken <- le.vn + 9L
        charleft <- floor(getOption("width") - chartaken - 4) # to have room for ...\n 
        cat(showColDataStruct(object, name, charleft))
        cat("\n")
      }
    }
  } 
)


#' @name data.struct
#' @description A constructor for data.struct from a dataframe
#' @export
data.struct <- function(...) {
  df <- data.frame(...)
  new("data.struct", ptr = DataFrameToDataStructR_(df), type = 0L)
}
