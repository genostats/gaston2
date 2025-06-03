a <- readBedFileMemory_("inst/extdata/LCT.bed", 503, 607)
A <- SNPMatrixToNumericMatrix(a)

# test GRM

g1 <- grm(a)

g2 <- tcrossprod(A) / ncol(A)
stopifnot( max(abs(g1 - g2)) < 5e-5 ) # ces calculs sont en float !

data("LCT", package = "gaston")
x <- gaston::as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

g3 <- gaston::GRM(x)*606/607
stopifnot( max(abs(g1 - g3)) < 5e-5 ) 


if(TRUE) {
RcppParallel::setThreadOptions(1)
set_num_thread(1)
mb1 <- microbenchmark::microbenchmark(grm(a), tcrossprod(A), gaston::GRM(x))

RcppParallel::setThreadOptions(4)
set_num_thread(4)
mb4 <- microbenchmark::microbenchmark(grm(a), gaston::GRM(x))

RcppParallel::setThreadOptions(8)
set_num_thread(8)
mb8 <- microbenchmark::microbenchmark(grm(a), gaston::GRM(x))
}

