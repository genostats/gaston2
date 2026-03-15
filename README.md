# Gaston 2

Gaston 2 is under development. The first release is planned for mid-May
2026.

## A short overview

The R package `gaston2` is a complete rewrite of `gaston`
<https://cran.r-project.org/package=gaston>.  
`Gaston` was first released in 2015 and works well for datasets of up to
10,000 individuals. `Gaston2` will allow users to manipulate larger
datasets without loading them into memory, by using memory-mapped
objects. It will also include additional features for dosage data.

The main features of `gaston2` will be:

  - Data manipulation and filtering for quality control.
  - Computation of genetic relationship matrices and principal component
    analysis (PCA).
  - Association studies (LMM and GLMM).

## A simple example

Here is a simple example of data manipulation with `gaston2`.

``` r
require(gaston2)
filename <- system.file("extdata", "LCT.bed", package="gaston2")

# reading in memory
x1 <- read.snp.matrix(filename)
x1
```

    ## A snp.matrix with 503 individuals and 607 SNPs
    ## Loaded in memory
    ## File: /home/rv/R/x86_64-pc-linux-gnu-library/4.5/gaston2/extdata/LCT.bed

``` r
# memory mapped object
x2 <- read.snp.matrix(filename, memory = FALSE)
x2
```

    ## A snp.matrix with 503 individuals and 607 SNPs
    ## On disk
    ## File: /home/rv/R/x86_64-pc-linux-gnu-library/4.5/gaston2/extdata/LCT.bed

``` r
# You can control wether the submatrices are in temporary files, or in memmory
x2[ 1:20, ]
```

    ## A snp.matrix with 20 individuals and 607 SNPs
    ## On disk
    ## File: /tmp/RtmpBsF4B6/gaston2e071b2b5931c7

``` r
x2[ 1:20, , "memory"]
```

    ## A snp.matrix with 20 individuals and 607 SNPs
    ## Loaded in memory

``` r
# Computing a LD matrix as an R object
A1 <- LD(x2, c(1, ncol(x2)), measure = "r2")
str(A1)
```

    ##  num [1:607, 1:607] 1 0.82913 0.08931 0.00554 0.43606 ...

``` r
# ...or as a memory mapped matrix, thanks to our package 'houba'
A2 <- LD(x2, c(1, ncol(x2)), measure = "r2", filename = tempfile())
A2
```

    ## A mmatrix with 607 rows and 607 cols
    ## data type:  float 
    ## Location: file  /tmp/RtmpBsF4B6/filee071b2f65d673 
    ## --- excerpt
    ##             [,1]       [,2]         [,3]         [,4]       [,5]
    ## [1,] 0.999996424 0.82913041 0.0893110931 0.0055414606 0.43606403
    ## [2,] 0.829130411 1.00000072 0.0742693171 0.0117575107 0.51423913
    ## [3,] 0.089311093 0.07426932 1.0000073910 0.0006292131 0.03957612
    ## [4,] 0.005541461 0.01175751 0.0006292131 1.0000023842 0.28635943
    ## [5,] 0.436064035 0.51423913 0.0395761244 0.2863594294 1.00000668
