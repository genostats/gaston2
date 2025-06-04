#include "SNPmatrix.h"
#include "extractIndsfromSNPmatrixMemory.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix> extractIndsfromSNPmatrixMemory_(Rcpp::XPtr<SNPmatrix> other, Rcpp::IntegerVector keep) {
  Rcpp::XPtr<SNPmatrix> pnewMat(new SNPmatrix);
  extractIndsfromSNPmatrixMemory(*other, keep, *pnewMat);
  return pnewMat;
}
