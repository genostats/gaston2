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

// [[Rcpp::export]]
Rcpp::DataFrame test_ds2() {
  std::vector<int> v;
  v.push_back(12);
  Column col(v);
  col.push_back(1);

  Column col2(datatype::DOUBLE);
  col2.push_back(3.14);
  col2.push_back(2.71);

  DataStruct DS;
  DS.push_back(col,  "a");
  DS.push_back(col2, "b");
  return DataStructToDataFrame(DS);
}
