#include <Rcpp.h>
#include <iostream>
#include "DataStruct.h"
#include "getDataStructPtr.h"

// [[Rcpp::export]]
Rcpp::IntegerVector dimDataStruct(Rcpp::S4 x) {
  DataStruct * pDS = getDataStructPtr(x);
  int nb_row = pDS->nrow();
  int nb_col = pDS->size();
  return Rcpp::IntegerVector::create(nb_row, nb_col);
}
