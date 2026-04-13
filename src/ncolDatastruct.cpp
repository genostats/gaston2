#include "DataStruct.h"
#include <iostream>
#include <Rcpp.h>


// [[Rcpp::export]]
size_t ncolDataStruct(Rcpp::XPtr<DataStruct> pDS) {
  return pDS->ncol();
}

