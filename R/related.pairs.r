#' @export
related.pairs <- function(K, min = 0.025, max = Inf) {
  n <- nrow(K)
  if(n != ncol(K)) stop("Not a square matrix")
  # K should be symmetric but we don't test that

  # extract indices of pairs between min and max
  # (hopefully) efficient way, that works both for matrices / mmatrices
  f <- function(x) which(x >= min & x <= max)
  L <- apply(K, 1, f, simplify = FALSE)
  L <- sapply(1:n, function(i) setdiff(L[[i]], 1:i))
  i <- rep(1:n, sapply(L, length))
  j <- unlist(L)

  D <- data.frame(i = i, j = j)

  # adding id_i id_j after in case user provided a matrix without dimnames
  D$id_i <- rownames(K)[i]
  D$id_j <- colnames(K)[j]
  
  # adding kinship coeff  
  D$k <- as.vector( K[i + (j - 1L)*n] )
  
  D 
}
