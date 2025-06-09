#include "SNPmatrix.h"
#include "SNPdosage.h"
#include "extractIndsfromDosagematrixMemory.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosage>> extractIndsfromDosagematrixMemory_(Rcpp::XPtr<SNPmatrix<SNPdosage>> other, Rcpp::IntegerVector keep) {
  Rcpp::XPtr<SNPmatrix<SNPdosage>> pnewMat(new SNPmatrix<SNPdosage>);
  extractIndsfromDosagematrixMemory(*other, keep, *pnewMat);
  return pnewMat;
}
