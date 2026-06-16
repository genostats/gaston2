#include <Rcpp.h>
#include <iostream>
#include "Column.h"
#include "DataStruct.h"
#include "getDataStructPtr.h"

// [[Rcpp::export]]
Rcpp::XPtr<DataStruct> extractDataStruct(Rcpp::S4 x, Rcpp::IntegerVector lines, Rcpp::IntegerVector columns) {
  DataStruct * pDS = getDataStructPtr(x);
  DataStruct * newDS(new DataStruct(*pDS, lines, columns));
  return Rcpp::XPtr<DataStruct>(newDS);
}
