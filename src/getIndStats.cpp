#include "SNPmatrix.h"
#include "SNPdosage.h"
#include "DataStruct.h"
#include <iostream>
#include <Rcpp.h>


Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame getIndStats(Rcpp::XPtr<SNPmatrix<>> pM, bool compute = false) {
  pM->compute_indStats(compute); // to force recomputing, if never accessed before, even with no force will compute
  return DataStructToDataFrame( pM->getIndStats() );
}

// [[Rcpp::export]]
Rcpp::S4 getIndStats_DataStruct(Rcpp::XPtr<SNPmatrix<>> pM, bool compute = false) {
    pM->compute_indStats(compute);
    Rcpp::XPtr<DataStruct> ptr(new DataStruct(pM->getIndStats()));
    Rcpp::S4 ds("data.struct");
    ds.slot("ptr") = ptr;
    return ds;
}

// [[Rcpp::export]]
Rcpp::DataFrame getIndStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool compute = false) {
  pM->compute_indStats(compute); // to force recomputing, if never accessed before, even with no force will compute
  return DataStructToDataFrame( pM->getIndStats() );
}


// [[Rcpp::export]]
Rcpp::S4 getIndStatsDosage_DataStruct(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool compute = false) {
    pM->compute_indStats(compute);
    Rcpp::XPtr<DataStruct> ptr(new DataStruct(pM->getIndStats()));
    Rcpp::S4 ds("data.struct");
    ds.slot("ptr") = ptr;
    return ds;
}

// [[Rcpp::export]]
void test_force_compute_indStats(Rcpp::XPtr<SNPmatrix<>> pM) {
  pM->compute_indStats(true);
}
