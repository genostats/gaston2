#' @useDynLib gaston2, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

.onAttach <- function(libname, pkgname) {
  rnorm(1); # force RNG seed initialisation (not done when the RNG is called from C++ code)

  n <- get_num_procs()
  n <- 2**max(0,floor(min(3, log2(n) - 1)))
  if(r_check_limit_cores())
    n <- 1
  set_num_threads(n)
  # packageStartupMessage("Gaston 2 set number of threads to ", n)
}

.onLoad <- function(libname, pkgname) {
  rnorm(1); # force RNG seed initialisation (not done when the RNG is called from C++ code)

  if(r_check_limit_cores())
    setThreadOptions(1)
}


# note : Il se pourrait qu'au lieu de renvoyer "true" certains systèmes renvoient "true " ou quelque chose comme ça
# le plus simple semble etre de tester que c'est défini (à autre chose que "false" au cas où, admettons)...
r_check_limit_cores <- function() { 
  Rcheck <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  (nchar(Rcheck[1]) > 0) & (Rcheck != "false")
}

