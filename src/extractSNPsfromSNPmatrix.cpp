#include "SNPmatrix.h"
#include <Rcpp.h>

//R exported function to have a matrix of only SNPs from "other" specified in "keep"
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix> extractSNPsfromSNPmatrix_(Rcpp::XPtr<SNPmatrix> other, Rcpp::IntegerVector keep) {
  Rcpp::XPtr<SNPmatrix> pnewMat(new SNPmatrix(*other, keep));
  return pnewMat;
}

