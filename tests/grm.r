library(gaston2)


a <- read.snp.matrix("inst/extdata/LCT.bed")

set.mode(a, "standardized.p")
A <- as.matrix(a)

# test GRM
g1 <- gaston2:::GRM(a)

g2 <- tcrossprod(A) / ncol(A)
stopifnot( max(abs(g1 - g2)) < 5e-5 ) # ces calculs sont en float !

data("LCT", package = "gaston")
x <- gaston::as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

g3 <- gaston::GRM(x)*606/607
stopifnot( max(abs(g1 - g3)) < 5e-5 ) 

# which(g1 != g3, arr.ind = TRUE)
# diff_res <- which( abs(g1 - g3) > 5e-5, arr.ind = TRUE) 
# equ_res <- which( abs(g1 - g3) < 5e-5, arr.ind = TRUE)
# b <- a[1:5, 1:5]
# g3 <- gaston::GRM(x[1:5,1:5])*606/607
# bg <- as.matrix(GRM(b, "test_grm"))

fn <- "grm_test"
g4 <- as.matrix(gaston2::GRM(a, filename=fn))
stopifnot( max(abs(g4 - g3)) < 5e-5 )

# tests with a calcul in floats
fn_float <- "grm_test_float"
g5 <- gaston2::GRM(a, filename=fn_float, usefloat=TRUE)
stopifnot( max(abs(as.matrix(g5) - g3)) < 5e-5 )

if(TRUE) {
RcppParallel::setThreadOptions(1)
gaston2:::set_num_thread(1)
mb1 <- microbenchmark::microbenchmark(gaston2:::GRM(a), gaston2:::GRM(a, filename=fn), gaston2:::GRM(a, filename=fn_float, usefloat=TRUE), tcrossprod(A), gaston::GRM(x))

RcppParallel::setThreadOptions(4)
gaston2:::set_num_thread(4)
mb4 <- microbenchmark::microbenchmark(gaston2:::GRM(a), gaston2:::GRM(a, filename=fn), gaston2:::GRM(a, filename=fn_float, usefloat=TRUE), gaston::GRM(x))

RcppParallel::setThreadOptions(8)
gaston2:::set_num_thread(8)
mb8 <- microbenchmark::microbenchmark(gaston2:::GRM(a), gaston2:::GRM(a, filename=fn), gaston2:::GRM(a, filename=fn_float, usefloat=TRUE), gaston::GRM(x))
}

#Cleanup grm_test
if (file.exists(fn)) {
  shut_up <- file.remove(fn)
}

#Cleanup grm_test_float
if (file.exists(fn_float)) {
  shut_up <- file.remove(fn_float)
}




