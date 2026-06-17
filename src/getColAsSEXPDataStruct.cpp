#include <Rcpp.h>
#include <iostream>
#include "RcppDataStruct.h"

// [[Rcpp::export]]
SEXP getColAsSEXPDataStruct(Rcpp::S4 x, std::string colname) {
  DataStruct * pDS = getDataStructPtr(x);
  Column & col = pDS->getColumn(colname);
  return ColumnToSEXP(col);
}

