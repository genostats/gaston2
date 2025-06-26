#' @name bind
#' @rdname bind
#' @title Binding methods for snp.matrix objects
#'
#' @aliases rbind2,snp.matrix
#' @aliases cbind2,snp.matrix
#'
#' @usage \S4method{rbind2}{snp.matrix,snp.matrix}(x,y,...)
#' @usage \S4method{cbind2}{snp.matrix,snp.matrix}(x,y,...)
#' TODO : check if missing version necessary ? 
#' @usage \S4method{cbind2}{snp.matrix,missing}(x,y,...)
#'
#' @param x,y snp.matrix objects
#' @param ... extra parameters (ignored)
#'
#' @description Methods allowing to use `cbind` and `rbind` with snp.matrix objects.
#'
#' @return A (new) snp.matrix matrix combining the arguments.
#'
#' @examples 
#' x <- dual( c(1, 3) )
#' y <- cbind(x, 2*x+1, 3*x+2, c(0,1))
#' y
#' d(y, "x1")
#' 
#'

# rbind
rbind2_snp_memory <- function(x, y, ...) {
  # print("Calling bindIndstoSNPmatrixMemory_")
  new_ptr <- bindIndstoSNPmatrixMemory_(x@ptr, y@ptr)
  new("snp.matrix", ptr = new_ptr, file = NULL, type = "memory")
}

#' @export
rbind2_snp_disk <- function(x, y, file, ...) {
  # add a check that file correct ?or let c++ do it ?
  # print("Calling bindIndstoSNPmatrixDisk_")
  new_ptr <- bindIndstoSNPmatrixDisk_(x@ptr, y@ptr, file)
  new("snp.matrix", ptr = new_ptr, file = file, type = "disk")
}


#' @export
setMethod("rbind2", c(x = "snp.matrix", y = "snp.matrix"),
  function(x, y, ..., file = NULL) {
    #print("using my rbind")
    if(...length() > 0)
      type <- match.arg(..1, c("disk", "memory"))
    else
      type <- x@type
      # sinon je prends le type de la première

    #get file from ... if it was passed ? 
    #else do an error
    # also modify the path with all prefix from r
    if (type == "disk") {
      if(type == "disk") {
        if(...length() > 1) 
          file <- ..2
        else
          file <- tempfile("gaston2")
      } else {
        file <- NULL
      }
    }
    if(type == "disk")
      rbind2_snp_disk(x, y, file = file, ...)
    else
      rbind2_snp_memory(x, y, ...)
  }
)

# for cbind,
# because the new matrix only exists as an "empty" object
# and all its content is simply pointed to 
# and stored in other matrix effectively handling it
# the type is a bit weird. 
# Could be defined as in memory without having any SNPs loaded
# but seems the most logical thing to do ?
# or else add a way to keep track of all files
# that it's refering to but seems a bit tedious...

#' @export
setMethod("cbind2", c(x = "snp.matrix", y = "snp.matrix"), 
  function(x, y, ...) {
  #print("Calling cbind_SNPmatrix")
  new_ptr <- cbind_SNPmatrix(x@ptr, y@ptr)
  new("snp.matrix", ptr = new_ptr, type = "memory")
  }
)