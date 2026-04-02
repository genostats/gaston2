#include "SNPmatrix.h"
#include "SNPdosage.h"
#include <iostream>
#include <Rcpp.h>


Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStats(Rcpp::XPtr<SNPmatrix<>> pM, bool compute = false) {
  // If never populated before (which should not happen), 
  // will create N0s, N1s... AND compute all SNPStats
  // else will return (if compute = false)
  pM->exportSNPStats(compute);
  return DataStructToDataFrame( pM->getSNPStats() );
}

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool compute = false) {
  // If never populated before (which should not happen), 
  // will create N0s, N1s... AND compute all SNPStats
  pM->exportSNPStats(compute);
  return DataStructToDataFrame( pM->getSNPStats() );
}
