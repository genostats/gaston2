#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>


Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStats(Rcpp::XPtr<SNPmatrix> pM) {
  return DataStructToDataFrame( pM->getSNPStats() );
}

