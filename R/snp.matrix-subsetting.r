#' @name extract
#' @title Extract parts of SNP matrix
#'
#' @description Extract columns / lines of a snp.matrix object
#'
#' @rdname extract
setMethod("[", c(x = "snp.matrix", i = "missing", j = "numeric", drop = "ANY"),
  function(x, i, j, ..., drop) {
    new_ptr <- extractSNPsfromSNPmatrix_(x@ptr, as.integer(j) - 1L)
    new("snp.matrix", ptr = new_ptr, file = NULL, type = x@type)
  }
)

#' @rdname extract
setMethod("[", c(x = "snp.matrix", i = "missing", j = "logical", drop = "ANY"),
  function(x, i, j, ..., drop) { 
    if(length(j) < ncol(x)) { # recycling
      j <- rep_len(j, ncol(x))
    }
    x[, which(j)]
  }
)

#' @rdname extract
setMethod("[", c(x = "snp.matrix", i = "numeric", j = "missing", drop = "ANY"),
  function(x, i, j, ..., drop) {
    if(...length() > 0)
      type <- match.arg(..1, c("disk", "memory"))
    else
      type <- x@type
   
    if(type == "disk") {
      if(...length() > 1) 
        file <- ..2
      else
        file <- tempfile("gaston2")
    } else {
      file <- NULL
    }

    if(type == "disk") 
      new_ptr <- extractIndsfromSNPmatrixDisk_(x@ptr, as.integer(i) - 1L, file)
    else
      new_ptr <- extractIndsfromSNPmatrixMemory_(x@ptr, as.integer(i) - 1L)

    new("snp.matrix", ptr = new_ptr, file = file, type = type)
  }
)

#' @rdname extract
setMethod("[", c(x = "snp.matrix", i = "logical", j = "missing", drop = "ANY"),
  function(x, i, j, ..., drop) {
    if(length(i) < nrow(x)) { # recycling
      i <- rep_len(i, nrow(x))
    }
    x[which(i), ]
  }
)

#' @rdname extract
setMethod("[", c(x = "snp.matrix", i = "index", j = "index", drop = "ANY"),
  function(x, i, j, ..., drop) {
    # on extrait d'abord les colonnes, moins couteux
    x[,j][i,,...]
  }
)

