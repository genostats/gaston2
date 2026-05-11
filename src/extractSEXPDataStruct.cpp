#include <Rcpp.h>

#include <iostream>

#include "Column.h"
#include "DataStruct.h"
#include "datatype.h"
#include <stdexcept> // for runtime_error

// [[Rcpp::export]]
SEXP extractSEXPDataStruct(Rcpp::XPtr<DataStruct> pDS, std::string colname) {
  if (!pDS) throw std::runtime_error("data.struct has a broken external pointer !");
  Column & col = pDS->getColumn(colname);
  datatype dt = col.type();
  if(dt == INT) 
    return Rcpp::wrap(*col.get<int>());
   else if(dt == DOUBLE)
    return Rcpp::wrap(*col.get<double>());
  else if(dt == FLOAT)
    return Rcpp::wrap(*col.get<float>());
  else if(dt == STRING)
    return Rcpp::wrap(*col.get<std::string>());
  else if(dt == BOOL)
    return Rcpp::wrap(*col.get<bool>());
  else
    Rcpp::stop("Unknown column type");
}
