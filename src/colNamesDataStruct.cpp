#include "DataStruct.h"
#include "getDataStructPtr.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
std::vector<std::string> colNamesDataStruct(Rcpp::S4 x) {
  DataStruct * pDS = getDataStructPtr(x);
  std::vector<std::string> names;
  size_t ncol = pDS->ncol();
  for(size_t i = 0; i < ncol; i++)
      names.push_back(pDS->colName(i));
  return names;
}
