#include "SNPmatrix.h"
#include "SNPdosage.h"
#include <iostream>
#include <Rcpp.h>


Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStats(Rcpp::XPtr<SNPmatrix<>> pM) {
  return DataStructToDataFrame( pM->getSNPStats() );
}

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM) {
  return DataStructToDataFrame( pM->getSNPStats() );
}
