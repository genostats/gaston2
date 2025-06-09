#' @name extract
#' @title Extract parts of DOSE matrix
#'
#' @description Extract columns / lines of a dose.matrix object
#'
#' @rdname extract
setMethod("[", c(x = "dose.matrix", i = "missing", j = "numeric", drop = "ANY"),
  function(x, i, j, ..., drop) {
    new_ptr <- extractSNPsfromDosagematrix_(x@ptr, as.integer(j) - 1L)
    new("dose.matrix", ptr = new_ptr, file = NULL, type = x@type)
  }
)

#' @rdname extract
setMethod("[", c(x = "dose.matrix", i = "missing", j = "logical", drop = "ANY"),
  function(x, i, j, ..., drop) { 
    if(length(j) < ncol(x)) { # recycling
      j <- rep_len(j, ncol(x))
    }
    x[, which(j)]
  }
)

#' @rdname extract
setMethod("[", c(x = "dose.matrix", i = "numeric", j = "missing", drop = "ANY"),
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
      new_ptr <- extractIndsfromDosagematrixDisk_(x@ptr, as.integer(i) - 1L, file)
    else
      new_ptr <- extractIndsfromDosagematrixMemory_(x@ptr, as.integer(i) - 1L)

    new("dose.matrix", ptr = new_ptr, file = file, type = type)
  }
)

#' @rdname extract
setMethod("[", c(x = "dose.matrix", i = "logical", j = "missing", drop = "ANY"),
  function(x, i, j, ..., drop) {
    if(length(i) < nrow(x)) { # recycling
      i <- rep_len(i, nrow(x))
    }
    x[which(i), ]
  }
)

#' @rdname extract
setMethod("[", c(x = "dose.matrix", i = "index", j = "index", drop = "ANY"),
  function(x, i, j, ..., drop) {
    # on extrait d'abord les colonnes, moins couteux
    x[,j][i,,...]
  }
)

