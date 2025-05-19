#include "Datastruct.h"
#include <iostream>
#include <Rcpp.h>
#include "debug.h"

DataStruct DataFrameToDataStruct(Rcpp::DataFrame DF);
Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);

// [[Rcpp::export]]
Rcpp::DataFrame test_ds(Rcpp::DataFrame DF, Rcpp::IntegerVector In) {
  DataStruct DS = DataFrameToDataStruct(DF);
  DataStruct DS2(DS, In);
  return DataStructToDataFrame(DS2);
}
