#' @export
ld.square <- function(snp.matrix, i1, i2, basename, usefloat = TRUE) {
  if(is.null(snp.matrix@ptr)) stop("Your snp.matrix has a broken ptr!")
  return(LD_square_mmatrix(snp.matrix@ptr, i1, i2, basename, TRUE, usefloat))
}

#' @export
ld.chunk <- function(snp.matrix, i1, i2, j1, j2, basename, usefloat = TRUE) {
  if(is.null(snp.matrix@ptr)) stop("Your snp.matrix has a broken ptr!")
  return(LD_chunk_mmatrix(snp.matrix@ptr, i1, i2, j1, j2, basename, TRUE, usefloat))
}

#' @export
ld.square.em <- function(snp.matrix, i1, i2, basename, usefloat = TRUE) {
  if(is.null(snp.matrix@ptr)) stop("Your snp.matrix has a broken ptr!")
  return(LD_square_EM_mmatrix(snp.matrix@ptr, i1, i2, basename, TRUE, usefloat))
}

#' @export
ld.chunk.em <- function(snp.matrix, i1, i2, j1, j2, basename, usefloat = TRUE) {
  if(is.null(snp.matrix@ptr)) stop("Your snp.matrix has a broken ptr!")
  return(LD_chunk_EM_mmatrix(snp.matrix@ptr, i1, i2, j1, j2, basename, TRUE, usefloat))
}
