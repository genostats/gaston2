require(gaston2)

a <- read.snp.matrix("inst/extdata/LCT.bed")
b <- read.snp.matrix("inst/extdata/LCT.bed", FALSE)

data("LCT", package = "gaston")
if(FALSE) {
# converting to integer matrix
A <- SNPMatrixToIntegerMatrix(a)

# checking values
cp <- table(A , LCT.gen, useNA = "always")
stopifnot( all(cp == diag(c(30286L, 69726L, 205306L, 3L))) )
}

# gaston bed matrix
x <- gaston::as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

# checking stat values 
cn <- c("N0", "N1", "N2", "NAs")
stopifnot( all(ind.stats(a)[,cn] == x@ped[,cn]) )

if(FALSE) {
# checking read fam file
a <- readBedFileMemory_("inst/extdata/LCT.bed", 503, 607)
readFamFile(a, "inst/extdata/LCT.fam")
head(getIndStats(a, FALSE))
head(getIndStats(a, TRUE))

readBimFile(a, "inst/extdata/LCT.bim")
head(getSNPStats(a))
}
