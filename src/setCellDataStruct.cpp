#include <Rcpp.h>
#include <iostream>
#include "RcppDataStruct.h"

// [[Rcpp::export]]
void setCellDataStruct(Rcpp::S4 x, int i, int j, SEXP value) {
  DataStruct * pDS = getDataStructPtr(x);
  Column & col = pDS->at(j);

  datatype dt = col.type();
  if (dt == INT)
    col.set<int>(i, Rf_asInteger(value));
  else if (dt == DOUBLE)
    col.set<double>(i, Rf_asReal(value));
  else if (dt == FLOAT)
    col.set<float>(i, (float) Rf_asReal(value));
  else if (dt == STRING)
    col.set<std::string>(i, Rcpp::as<std::string>(value));
  else if (dt == BOOL)
    col.set<bool>(i, Rf_asBool(value));
  else
    Rcpp::stop("Unknown column type");
}

