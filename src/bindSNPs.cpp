#include "SNPmatrix.h"
#include "SNPdosage.h"
#include <Rcpp.h>

//R exported function to have a new SNPmatrix<SNPvector> of combined snps from
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<>> cbind_SNPmatrix(Rcpp::XPtr<SNPmatrix<>> first_matrix, Rcpp::XPtr<SNPmatrix<>> scd_matrix) {
  Rcpp::XPtr<SNPmatrix<>> pnewMat(new SNPmatrix<>(*first_matrix, *scd_matrix));
  return pnewMat;
}

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosage>> cbind_Dosagematrix(Rcpp::XPtr<SNPmatrix<SNPdosage>> first_matrix, Rcpp::XPtr<SNPmatrix<SNPdosage>> scd_matrix) {
  Rcpp::XPtr<SNPmatrix<SNPdosage>> pnewMat(new SNPmatrix<SNPdosage>(*first_matrix, *scd_matrix));
  return pnewMat;
}