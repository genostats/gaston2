#include "SNPmatrix.h"
#include "SNPdosage.h"
#include "extractIndsfromDosagematrixDisk.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosage>> extractIndsfromDosagematrixDisk_(Rcpp::XPtr<SNPmatrix<SNPdosage>> other, Rcpp::IntegerVector keep, std::string path_str) {
  Rcpp::XPtr<SNPmatrix<SNPdosage>> pnewMat(new SNPmatrix<SNPdosage>);
  extractIndsfromDosagematrixDisk(*other, keep, path_str, *pnewMat);
  return pnewMat;
}
