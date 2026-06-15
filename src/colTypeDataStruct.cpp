#include "DataStruct.h"
#include "getDataStructPtr.h"
#include "Column.h"
#include "datatype.h" // for typeToString
#include <iostream>
#include <stdexcept>// for runtime_error exception
#include <Rcpp.h>

// [[Rcpp::export]]
std::string colTypeDataStruct(Rcpp::S4 x, std::string colName) {
  DataStruct * pDS = getDataStructPtr(x);
  std::string type;
  if (pDS->hasColumn(colName)){
    // copie pas génante parce que uniquement le ptr, pas vect de datas
    Column col = pDS->getColumn(colName);
    type = typeToString(col.type());
  } else {
    Rcpp::stop("Column doesn't exist");
  }
  return type;
}
