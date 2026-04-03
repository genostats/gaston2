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

all(c("N0", "N1", "N2", "NAs") %in% colnames(snp.stats(a_ref)))
all(snp.stats(a_ref) == snp.stats(a_test))
all(c("N0", "N1", "N2", "NAs") %in% colnames(ind.stats(a_ref)))
all(ind.stats(a_ref) == ind.stats(a_test))
all(as.matrix(a_ref) == as.matrix(a_test))

# TO CHANGE WITH YOUR OWN DOSAGE FILE !!
d <- read.dose.matrix("truc", FALSE)
d1 <- d[,26]
d2 <- d[,27]
d_test <- cbind(d1, d2)

d_ref <- d[,26:27]

all(c("N0", "N1", "N2", "NAs") %in% colnames(snp.stats(d_ref)))
all(snp.stats(d_ref) == snp.stats(d_test))
all(c("N0", "N1", "N2", "NAs") %in% colnames(ind.stats(d_ref)))
all(ind.stats(d_ref) == ind.stats(d_test))
all(as.matrix(d_ref) == as.matrix(d_test))




# -------- testing the rbind (append of inds)



rb1 <- a[1:100,]
rb2 <- a[101:200]
#rb_test <- rbind(rb1, rb2)
#i cannot pass argument to rbind, so I'm having this roundaround
rb_test <- rbind2_snp_disk(rb1, rb2, "tmp/test_rbind")
#line needed to compare with ref (or else not calling it on ref either)

rbind_ref <- a[1:200,]

all(c("N0", "N1", "N2", "NAs") %in% colnames(snp.stats(rbind_ref)))
all(snp.stats(rbind_ref) == snp.stats(rb_test))
all(c("N0", "N1", "N2", "NAs") %in% colnames(ind.stats(rbind_ref)))
all(ind.stats(rbind_ref) == ind.stats(rb_test))
all(as.matrix(rbind_ref) == as.matrix(rb_test), na.rm = TRUE)

# to update w/ ur dos file
d_disk <- read.dose.matrix("truc", FALSE)
d <- read.dose.matrix("truc")
d1 <- d[4:6,]
d2 <- d_disk[7:8,] #having a mismatch, combining on memory and on disk

d_test <- rbind(d1, d2)
disk_test <- rbind2_dose_disk(d1, d2, "tmp/test_rbind.dosf")

d_ref <- d[4:8,]

all(c("N0", "N1", "N2", "NAs") %in% colnames(snp.stats(d_ref)))
all(snp.stats(d_ref) == snp.stats(d_test))
all(snp.stats(d_ref) == snp.stats(disk_test))
all(c("N0", "N1", "N2", "NAs") %in% colnames(ind.stats(rbind_ref)))
all(ind.stats(d_ref) == ind.stats(d_test))
all(ind.stats(d_ref) == ind.stats(disk_test))
all(as.matrix(d_ref) == as.matrix(d_test))
all(as.matrix(d_ref) == as.matrix(disk_test))



# -------- testing a bit of both

e <- read.snp.matrix("inst/extdata/LCT.bed")
e1 <- e[2:4,1:8]
  
data("LCT", package = "gaston")
eref <- gaston::as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

e1ref <- eref[2:4,1:8]

all(c("N0", "N1", "N2", "NAs") %in% colnames(snp.stats(e1)))
all(snp.stats(e1) == e1ref@snps[c("chr", "id", "dist", "pos", "A1", "A2", "N0", "N1", "N2", "NAs")])
all(c("N0", "N1", "N2", "NAs") %in% colnames(ind.stats(e1)))
all(ind.stats(e1) == e1ref@ped[c("famid", "id", "father", "mother", "sex", "pheno", "N0", "N1", "N2", "NAs")])
# all(as.matrix(e1) == as.matrix(e1ref))

#eventually add a test on dosages, but no easy comparison with gaston

#Cleanup files leftover from test
if (file.exists("tmp/test_rbind")) {
  file.remove("tmp/test_rbind")
}

if (file.exists("tmp/test_rbind.dosf")) {
  file.remove("tmp/test_rbind.dosf")
}
