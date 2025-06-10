require(gaston2)

a <- read.snp.matrix("inst/extdata/LCT.bed")

w <- rep(TRUE, 607)
gaston2:::LD_thin_(a@ptr, 0.2, 1000, 0, w)


# gaston bed matrix
data("LCT", package = "gaston")
x <- gaston::as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

R <- gaston2:::LD_square(a@ptr, 0, 606)
w1 <- gaston::LD.thin(x, 0.2, 1000, extract = FALSE, keep = "le")

stopifnot(all(w == w1))


