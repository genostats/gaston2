require(gaston2)


# -------- testing the cbind (append of snps)
a <- read.snp.matrix("inst/extdata/LCT.bed")
a1 <- a[,1]
a2 <- a[,2]

a_test <- gaston2:::cbind_SNPmatrix(a1@ptr, a2@ptr)
mtx_a_test <-new("snp.matrix", ptr = a_test, file = NULL, type = "memory")
#snp.stats(mtx_a_test)
gaston2:::getIndStats(mtx_a_test@ptr, FALSE)
#ind.stats(matrix_a_test) #calls getIndStats with compute == true, so N0, N1...


a_ref <- a[,1:2]
gaston2:::exportSNPStats(a_ref@ptr)


all(snp.stats(a_ref) == snp.stats(mtx_a_test))
all(ind.stats(a_ref) == ind.stats(mtx_a_test))
all(as.matrix(a_ref) == as.matrix(mtx_a_test))

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




# -------- testing the rbind (append of snps)

rb1 <- a[1:100,]
rb2 <- a[101:200]

rb_test <- gaston2:::bindIndstoSNPmatrixDisk_(rb1@ptr, rb2@ptr, "tmp/test_rbind")
rb_mtx_test <-new("snp.matrix", ptr = rb_test, file = "tmp/test_rbind", type = "disk")
#snp.stats(rb_mtx_test)
#gaston2:::getIndStats(rb_mtx_test@ptr, FALSE)
#ind.stats(rb_mtx_test) #calls getIndStats with compute == true, so N0, N1...

#line needed to compare with ref (or else not calling it on ref either)
gaston2:::exportSNPStats(rb_mtx_test@ptr)

rbind_ref <- a[1:200,]
gaston2:::exportSNPStats(rbind_ref@ptr)


all(snp.stats(rbind_ref) == snp.stats(rb_mtx_test))
all(ind.stats(rbind_ref) == ind.stats(rb_mtx_test))
all(as.matrix(rbind_ref) == as.matrix(rb_mtx_test), na.rm = TRUE)

# to update w/ dos file
d_disk <- read.dose.matrix("truc", FALSE)
d <- read.dose.matrix("truc")
d1 <- d[1,]
d2 <- d_disk[2,] #having a mismatch, combining on memory and on disk

d_test <- gaston2:::bindIndstoDosagematrixMemory_(d1@ptr, d2@ptr)
d_test_disk <- gaston2:::bindIndstoDosagematrixDisk_(d1@ptr, d2@ptr, "tmp/test_rbind.dosf")

#same here
matrix_d_test <-new("dose.matrix", ptr = d_test, file = NULL, type = "memory")

matrix_disk_test <-new("dose.matrix", ptr = d_test_disk, file = "tmp/test_rbind.dosf", type = "disk")
snp.stats(matrix_d_test)
gaston2:::getIndStats(matrix_d_test@ptr, FALSE)
ind.stats(matrix_d_test)

gaston2:::exportSNPStats(matrix_d_test@ptr)
gaston2:::exportSNPStats(matrix_disk_test@ptr)

d_ref <- d[1:2,]
gaston2:::exportSNPStats(d_ref@ptr)

all(snp.stats(d_ref) == snp.stats(matrix_d_test))
all(snp.stats(d_ref) == snp.stats(matrix_disk_test))
all(ind.stats(d_ref) == ind.stats(matrix_d_test))
all(ind.stats(d_ref) == ind.stats(matrix_disk_test))
all(as.matrix(d_ref) == as.matrix(matrix_d_test))
all(as.matrix(d_ref) == as.matrix(matrix_disk_test))
# 
