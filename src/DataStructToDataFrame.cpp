#include "Datastruct.h"
#include <iostream>
#include <Rcpp.h>


SEXP ColumnToSEXP(const Column & col) {
  datatype dt = col.type();
  if(dt == INT) 
    return Rcpp::wrap(*col.get<int>());
   else if(dt == DOUBLE)
    return Rcpp::wrap(*col.get<double>());
  else if(dt == FLOAT)
    return Rcpp::wrap(*col.get<float>());
  else if(dt == STRING)
    return Rcpp::wrap(*col.get<std::string>()); 
  else
    Rcpp::stop("Unknown column type");
}

Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS) {
  Rcpp::DataFrame DF;
  for(unsigned int i = 0; i < DS.size(); i++) {
    DF.push_back( ColumnToSEXP(DS.at(i)), DS.colName(i) );
  }
  return DF;
}


