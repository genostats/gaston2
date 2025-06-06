#include "SNPmatrix.h"
#include "extractIndsfromSNPmatrixDisk.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<>> extractIndsfromSNPmatrixDisk_(Rcpp::XPtr<SNPmatrix<>> other, Rcpp::IntegerVector keep, std::string path_str) {
  Rcpp::XPtr<SNPmatrix<>> pnewMat(new SNPmatrix<>);
  extractIndsfromSNPmatrixDisk(*other, keep, path_str, *pnewMat);
  return pnewMat;
}
