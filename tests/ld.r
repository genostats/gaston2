require(snipsnop)

ctl <- c(0.999999999999998, -0.0422681970547771, 0.916897472541938, 
8.21480811571587e-05, 0.994618132256276, 0.917510904376256, -0.0422681970547771, 
1, -0.0495072232188823, -0.00207325180667775, -0.0384700353407459, 
-0.048758086954353, 0.916897472541938, -0.0495072232188823, 0.999999999999997, 
0.0149290087967859, 0.917568304879441, 0.999508789301337, 8.21480811571587e-05, 
-0.00207325180667775, 0.0149290087967859, 0.999999999999993, 
-0.000985385235842333, 0.0151495088836487, 0.994618132256276, 
-0.0384700353407459, 0.917568304879441, -0.000985385235842333, 
1.00000000000001, 0.918190151673646, 0.917510904376256, -0.048758086954353, 
0.999508789301337, 0.0151495088836487, 0.918190151673646, 0.999999999999994)


a <- readBedFileMemory_("inst/extdata/LCT.bed", 503, 607)
computeSNPStats(a)

# ----- checking LD_square ----
M <- LD_square(a, 167, 172)
stopifnot( max(abs(M-ctl)) < 1e-10 )

# checking the 5 cases of LD_chunk
M0 <- LD_chunk(a, 167, 170, 170, 172)
stopifnot( all(M[1:4, 4:6] == M0) )

M1 <- LD_chunk(a, 167, 171, 170, 172)
stopifnot( all(M[1:5, 4:6] == M1) )

M2 <- LD_chunk(a, 169, 172, 167, 171)
stopifnot( all(M[3:6, 1:5] == M2) )

M3 <- LD_chunk(a, 169, 171, 167, 172)
stopifnot( all(M[3:5, 1:6] == M3) )

M4 <- LD_chunk(a, 167, 172, 168, 171)
stopifnot( all(M[1:6, 2:5] == M4) )


# ---- checking EM version ----
D <- LD_square_EM(a, 167, 172) 

D0 <- LD_chunk_EM(a, 167, 170, 170, 172)
stopifnot( all(D[1:4, 4:6] == D0) )

D1 <- LD_chunk_EM(a, 167, 171, 170, 172)
stopifnot( all(D[1:5, 4:6] == D1) )

D2 <- LD_chunk_EM(a, 169, 172, 167, 171)
stopifnot( all(D[3:6, 1:5] == D2) )

D3 <- LD_chunk_EM(a, 169, 171, 167, 172)
stopifnot( all(D[3:5, 1:6] == D3) )

D4 <- LD_chunk_EM(a, 167, 172, 168, 171)
stopifnot( all(D[1:6, 2:5] == D4) )

# --- checking big memory

LD_square_bigmemory(a, 0, 606, "essai.bm")
x <- attach.big.matrix("essai.bm.desc")
y <- LD_square(a, 0, 606)
range(y - as.matrix(x))




if(FALSE) {
data("LCT", package = "gaston")
x <- gaston::as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

mb1 <- microbenchmark::microbenchmark(M1 <- .Call(gaston:::`_gaston_LD`, PACKAGE = "gaston", x@bed, x@mu, x@sigma, 0, 606), M2 <- LD_square(a, 0, 606), times = 40)
stopifnot( max(abs(M1- M2)) < 1e-12 )

mb2 <- microbenchmark::microbenchmark(D <- LD_square_EM(a, 0, 606), times = 40)
}
