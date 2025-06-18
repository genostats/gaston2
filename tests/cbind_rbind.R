require(gaston2)


# -------- testing the cbind (append of snps)
a <- read.snp.matrix("inst/extdata/LCT.bed")
a1 <- a[,1]
a2 <- a[,2]

a_test <- gaston2:::cbind_SNPmatrix(a1@ptr, a2@ptr)
matrix_a_test <-new("snp.matrix", ptr = a_test, file = NULL, type = "memory")
snp.stats(matrix_a_test)
gaston2:::getIndStats(matrix_a_test@ptr, FALSE)
ind.stats(matrix_a_test) #calls getIndStats with compute == true, so N0, N1...


a_ref <- a[,1:2]
gaston2:::exportSNPStats(a_ref@ptr)


all(snp.stats(a_ref) == snp.stats(matrix_a_test))
all(ind.stats(a_ref) == ind.stats(matrix_a_test))
all(as.matrix(a_ref) == as.matrix(matrix_a_test))


# TO CHANGE WITH YOUR OWN DOSAGE FILE !!
d <- read.dose.matrix("truc", FALSE)
d1 <- d[,1]
d2 <- d[,2]

d_test <- gaston2:::cbind_Dosagematrix(d1@ptr, d2@ptr)
#same here !
matrix_d_test <-new("dose.matrix", ptr = d_test, file = "truc.dosf", type = "disk")
snp.stats(matrix_d_test)
gaston2:::getIndStats(matrix_d_test@ptr, FALSE)
ind.stats(matrix_d_test)

d_ref <- d[,1:2]
gaston2:::exportSNPStats(d_ref@ptr)

all(snp.stats(d_ref) == snp.stats(matrix_d_test))
all(ind.stats(d_ref) == ind.stats(matrix_d_test))
all(as.matrix(d_ref) == as.matrix(matrix_d_test))

