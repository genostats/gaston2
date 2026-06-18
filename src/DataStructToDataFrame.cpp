#include <Rcpp.h>
#include <iostream>
#include "RcppDataStruct.h"

Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS, Rcpp::IntegerVector I) {
  Rcpp::DataFrame DF;
  for(int i : I) {
    DF.push_back( ColumnToSEXP(DS.at(i)), DS.colName(i) );
  }
  return DF;
}
