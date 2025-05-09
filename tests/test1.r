require(snipsnop)
require(gaston)
x <- readBedFileMemory("inst/extdata/LCT.bed", 503, 607)
A <- SNPMatrixToIntegerMatrix(x)
cp <- table(A , LCT.gen, useNA = "always")

stopifnot( all(cp == diag(c(30286L, 69726L, 205306L, 3L))) )
