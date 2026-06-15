#include <Rcpp.h>
#include <iostream>
#include "Column.h"
#include "DataStruct.h"
#include "getDataStructPtr.h"
#include "datatype.h"
#include <stdexcept> // for runtime_error

/* Extract specified Columns (with names) from a DataStruct 
and pushes them into a new one before returning. 
*/

// [[Rcpp::export]]
Rcpp::XPtr<DataStruct> extractColDataStruct_(Rcpp::S4 x, std::vector<std::string> colnames) {
  DataStruct * pDS = getDataStructPtr(x);

  DataStruct * newDS = new DataStruct();
  int sizenames = colnames.size();
  for(size_t i = 0; i < sizenames; i++) {
    newDS->push_back(pDS->getColumn(colnames[i]), colnames[i]);
  }
  return Rcpp::XPtr<DataStruct>(newDS); 
}
