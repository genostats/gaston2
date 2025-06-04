#' @name snp.matrix
#' @rdname snp.matrix
#' @title SNP matrix
#' 
#' @description A class for manipulation of genotype data
#' @exportClass snp.matrix
setClass("snp.matrix", slots = c(ptr = "externalptr", file = "characterOrNull", type = "character"))

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


