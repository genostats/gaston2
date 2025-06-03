#' @title SNP matrix
#' 
#' @exportClass snp.matrix
#'
setClass("snp.matrix", slots = c(ptr = "externalptr", file = "character", type = "character"))

setMethod("show", "snp.matrix", 
  function(object) {
    if(isnullptr(object@ptr)) {
      cat("A snp.matrix with a broken external ptr\n")
    } else {
      cat("A snp.matrix with", nrow(object), "individuals and", ncol(object), "SNPs\n")
      if(object@type == "memory") 
        cat("Loaded in memory\n")
      if(object@type == "disk")
        cat("On disk\n")
      if(!is.null(object@file)) 
        cat("File:", object@file, "\n")
    }
  } 
)


