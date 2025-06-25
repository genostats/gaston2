require(gaston2)


# -------- testing the cbind (append of snps)
a <- read.snp.matrix("inst/extdata/LCT.bed")
a1 <- a[,1]
a2 <- a[,2]
a_test <- cbind(a1, a2)
# internal calls to that :
# a_test <- gaston2:::cbind_SNPmatrix(a1@ptr, a2@ptr)
# mtx_a_test <-new("snp.matrix", ptr = a_test, file = NULL, type = "memory")

a_ref <- a[,1:2]
gaston2:::exportSNPStats(a_ref@ptr)


all(snp.stats(a_ref) == snp.stats(a_test))
all(ind.stats(a_ref) == ind.stats(a_test))
all(as.matrix(a_ref) == as.matrix(a_test))

# TO CHANGE WITH YOUR OWN DOSAGE FILE !!
d <- read.dose.matrix("truc", FALSE)
d1 <- d[,26]
d2 <- d[,27]
d_test <- cbind(d1, d2)

d_ref <- d[,26:27]
gaston2:::exportSNPStats(d_ref@ptr)

all(snp.stats(d_ref) == snp.stats(d_test))
all(ind.stats(d_ref) == ind.stats(d_test))
all(as.matrix(d_ref) == as.matrix(d_test))




# -------- testing the rbind (append of snps)



rb1 <- a[1:100,]
rb2 <- a[101:200]
#rb_test <- rbind(rb1, rb2)
#i cannot pass argument to rbind, so I'm having this roundaround
rb_test <- rbind2_snp_disk(rb1, rb2, "tmp/test_rbind")
#line needed to compare with ref (or else not calling it on ref either)
gaston2:::exportSNPStats(rb_test@ptr)

rbind_ref <- a[1:200,]
gaston2:::exportSNPStats(rbind_ref@ptr)


all(snp.stats(rbind_ref) == snp.stats(rb_test))
all(ind.stats(rbind_ref) == ind.stats(rb_test))
all(as.matrix(rbind_ref) == as.matrix(rb_test), na.rm = TRUE)

# to update w/ ur dos file
d_disk <- read.dose.matrix("truc", FALSE)
d <- read.dose.matrix("truc")
d1 <- d[4:6,]
d2 <- d_disk[7:8,] #having a mismatch, combining on memory and on disk

d_test <- rbind(d1, d2)
disk_test <- rbind2_dose_disk(d1, d2, "tmp/test_rbind.dosf")

gaston2:::exportSNPStats(d_test@ptr)
gaston2:::exportSNPStats(disk_test@ptr)

d_ref <- d[4:8,]
gaston2:::exportSNPStats(d_ref@ptr)

all(snp.stats(d_ref) == snp.stats(d_test))
all(snp.stats(d_ref) == snp.stats(disk_test))
all(ind.stats(d_ref) == ind.stats(d_test))
all(ind.stats(d_ref) == ind.stats(disk_test))
all(as.matrix(d_ref) == as.matrix(d_test))
all(as.matrix(d_ref) == as.matrix(disk_test))