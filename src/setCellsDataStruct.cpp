#include <Rcpp.h>
#include <iostream>
#include "RcppDataStruct.h"

// [[Rcpp::export]]
void setCellsDataStruct(Rcpp::S4 x, Rcpp::IntegerVector I, int j, SEXP value) {
  DataStruct * pDS = getDataStructPtr(x);
  Column & col = pDS->at(j);

  size_t n = I.size();
  datatype dt = col.type();
  if (dt == INT) {
    Rcpp::IntegerVector Rvec(value);
    size_t m = Rvec.size();
    for(int k = 0; k < n; k++) col.set<int>(I[k], Rvec[k % m]);
  } else if (dt == DOUBLE) {
    Rcpp::NumericVector Rvec(value);
    size_t m = Rvec.size();
    for(int k = 0; k < n; k++) col.set<double>(I[k], Rvec[k % m]);
  } else if (dt == FLOAT) {
    Rcpp::NumericVector Rvec(value);
    size_t m = Rvec.size();
    for(int k = 0; k < n; k++) col.set<float>(I[k], (float) Rvec[k % m]);
  } else if (dt == STRING) {
    Rcpp::CharacterVector Rvec(value);
    size_t m = Rvec.size();
    for(int k = 0; k < n; k++) col.set<std::string>(I[k], Rcpp::as<std::string>(Rvec[k % m]));
  } else if (dt == BOOL) {
    Rcpp::LogicalVector Rvec(value);
    size_t m = Rvec.size();
    for(int k = 0; k < n; k++) col.set<bool>(I[k], (bool) Rvec[k % m]);
  } else {
    Rcpp::stop("Unknown column type");
  }
}

