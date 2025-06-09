#include "SNPmatrix.h"
#include "SNPdosage.h"
#include <Rcpp.h>

//R exported function to have a dosage matrix of only SNPs from "other" specified in "keep"
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosage>> extractSNPsfromDosagematrix_(Rcpp::XPtr<SNPmatrix<SNPdosage>> other, Rcpp::IntegerVector keep) {
  Rcpp::XPtr<SNPmatrix<SNPdosage>> pnewMat(new SNPmatrix<SNPdosage>(*other, keep));
  return pnewMat;
}
