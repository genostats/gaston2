#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>


Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame getIndStats(Rcpp::XPtr<SNPmatrix> pM) {
  pM->compute_indStats();
  return DataStructToDataFrame( pM->getIndStats() );
}

