#include <Rcpp.h>

#include <iostream>

#include "Column.h"
#include "DataStruct.h"
#include "datatype.h"
#include <stdexcept> // for runtime_error

/* Extract specified Columns (with names) from a DataStruct 
and pushes them into a new one before returning. 
*/

DataStruct *extractcolDataStruct(const DataStruct & DS, std::vector<std::string> colnames ) {
  DataStruct * newDS = new DataStruct();
  int sizenames = colnames.size();

  for(size_t i = 0; i < sizenames; i++) {
    newDS->push_back(DS.getColumn(colnames[i]), colnames[i]);
  }
  return newDS;
}

// [[Rcpp::export]]
Rcpp::XPtr<DataStruct> extractcolDataStruct_(Rcpp::XPtr<DataStruct> ogDS, std::vector<std::string> colnames ) {
  Rcpp::XPtr<DataStruct> pDS(extractcolDataStruct(*ogDS, colnames));
  return pDS; 
}
