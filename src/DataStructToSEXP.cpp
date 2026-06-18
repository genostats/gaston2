#include <Rcpp.h>
#include <iostream>
#include "RcppDataStruct.h"

// sends a DataFrame or a List
// [[Rcpp::export]]
SEXP DataStructToSEXP(Rcpp::S4 x, Rcpp::IntegerVector I) {
  DataStruct * pDS = getDataStructPtr(x);
  size_t nb_col = pDS->size();

  // empty DataFrame
  if(I.size() == 0) return Rcpp::DataFrame();

  // check of if all column have same size
  size_t ref_size = pDS->at(I(0)).size();
  bool all_equal = true;
  for(int i : I) {
    if (pDS->at(i).size() != ref_size) {
      all_equal = false;
      break;
    }
  }

  if(all_equal) // all same size, DF possible
    return DataStructToDataFrame(*pDS, I);
  else // fall back to list
    return DataStructToList(*pDS, I);
}
