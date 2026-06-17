#include <Rcpp.h>
#include <iostream>
#include "RcppDataStruct.h"

// sends a DataFrame or a List
// [[Rcpp::export]]
SEXP DataStructToSEXP(Rcpp::S4 x) {
  DataStruct * pDS = getDataStructPtr(x);
  size_t nb_col = pDS->size();
  // First need to run a check of if all column have same size
  size_t ref_size = pDS->at(0).size();
  bool all_equal = true;
  for(unsigned int i = 1; i < nb_col; i++) {
    if (pDS->at(i).size() != ref_size) {
      all_equal = false;
      break;
    }
  }
  if(all_equal) // all same size, DF possible
    return DataStructToDataFrame(*pDS);
  else // fall back to list
    return DataStructToList(*pDS);
}
