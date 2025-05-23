require(snipsnop)
require(gaston)
a <- readBedFileMemory_("inst/extdata/LCT.bed", 503, 607)
A <- SNPMatrixToIntegerMatrix(a)
cp <- table(A , LCT.gen, useNA = "always")

stopifnot( all(cp == diag(c(30286L, 69726L, 205306L, 3L))) )

x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

cn <- c("N0", "N1", "N2", "NAs")
stopifnot( all(getIndStats(a)[,cn] == x@ped[,cn]) )

