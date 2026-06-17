#include <Rcpp.h>
#include <iostream>
#include "DataStruct.h"
#include "getDataStructPtr.h"

// [[Rcpp::export]]
int nrowDataStruct(Rcpp::S4 x) {
  DataStruct * pDS = getDataStructPtr(x);
  return pDS->nrow();
}
