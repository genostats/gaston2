#include "SNPmatrix.h"
#include "SNPdosage.h"
#include <iostream>
#include <Rcpp.h>


Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame getIndStats(Rcpp::XPtr<SNPmatrix<>> pM, bool compute = false) {
  if (compute) pM->compute_indStats(compute); // to force recomputing, if never accessed before, even with no force will compute
  return DataStructToDataFrame( pM->getIndStats() );
}

// [[Rcpp::export]]
Rcpp::DataFrame getIndStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool compute = false) {
  if (compute) pM->compute_indStats(compute); // to force recomputing, if never accessed before, even with no force will compute
  return DataStructToDataFrame( pM->getIndStats() );
}

// [[Rcpp::export]]
void test_force_compute_indStats(Rcpp::XPtr<SNPmatrix<>> pM) {
  pM->compute_indStats(true);
}
