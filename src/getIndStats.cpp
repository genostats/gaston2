#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>


Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame getIndStats(Rcpp::XPtr<SNPmatrix> pM, bool compute = true) {
  if(compute) pM->compute_indStats();
  return DataStructToDataFrame( pM->getIndStats() );
}

// [[Rcpp::export]]
void test_force_compute_indStats(Rcpp::XPtr<SNPmatrix> pM) {
  pM->compute_indStats(true);
}
