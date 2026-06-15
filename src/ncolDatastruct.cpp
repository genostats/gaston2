#include "DataStruct.h"
#include "getDataStructPtr.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
size_t ncolDataStruct(Rcpp::S4 x) {
  DataStruct * pDS = getDataStructPtr(x);
  return pDS->ncol();
}

