#include <Rcpp.h>

#include <iostream>

#include "DataStruct.h"
#include "SNPdosage.h"
#include "SNPmatrix.h"

Rcpp::DataFrame DataStructToDataFrame(const DataStruct& DS);

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStats(Rcpp::XPtr<SNPmatrix<>> pM, bool compute = false) {
  // If never exported before will create N0s, N1s... AND compute all SNPStats
  // else will return (if compute = false)
  pM->exportSNPStats(compute);
  return DataStructToDataFrame(pM->getSNPStats());
}

// [[Rcpp::export]]
Rcpp::S4 getSNPStats_DataStruct(Rcpp::XPtr<SNPmatrix<>> pM, bool compute = false) {
  pM->exportSNPStats(compute);
  Rcpp::XPtr<DataStruct> ptr(new DataStruct(pM->getSNPStats()));
  Rcpp::S4 ds("data.struct");
  ds.slot("ptr") = ptr;
  return ds;
}

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool compute = false) {
  // If never exported before will create N0s, N1s... AND compute all SNPStats
  pM->exportSNPStats(compute);
  return DataStructToDataFrame(pM->getSNPStats());
}

// [[Rcpp::export]]
Rcpp::DataFrame getSNPStatsDosage_DataStruct(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool compute = false) {
  // If never exported before will create N0s, N1s... AND compute all SNPStats
  pM->exportSNPStats(compute);
  Rcpp::XPtr<DataStruct> ptr(new DataStruct(pM->getSNPStats()));
  Rcpp::S4 ds("data.struct");
  ds.slot("ptr") = ptr;
  return ds;
}
