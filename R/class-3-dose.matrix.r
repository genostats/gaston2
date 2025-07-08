#' @name dose.matrix
#' @rdname dose.matrix
#' @title Dosage matrix
#' 
#' @description A class for manipulation of genotype dosages
#' @exportClass dose.matrix
setClass("dose.matrix", slots = c(ptr = "externalptr", file = "characterOrNull", type = "character"))

setMethod("show", "dose.matrix", 
  function(object) {
    if(isnullptr(object@ptr)) {
      cat("A dose.matrix with a broken external ptr\n")
    } else {
      cat("A dose.matrix with", nrow(object), "individuals and", ncol(object), "SNPs\n")
      if(object@type == "memory") 
        cat("Loaded in memory\n")
      if(object@type == "disk")
        cat("On disk\n")
      if(!is.null(object@file)) 
        cat("File:", object@file, "\n")
    }
  } 
)

# de ce que je comprends pas besoin du setGeneric again, compilé à class-2-snp.matrix
#' @export
setMethod("add.descriptor.file", "dose.matrix", 
function(object) {
  # no need to check if memory or disk ?
  # TODO : add check NULL on file
  # RV's function : mk.descriptor.file <- function(path, nrow, ncol, type)
  mk.descriptor.file(object@file, nrow(object), ncol(object), "double")
})