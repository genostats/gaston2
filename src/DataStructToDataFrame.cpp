#include <Rcpp.h>
#include <iostream>
#include "RcppDataStruct.h"

Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS) {
  Rcpp::DataFrame DF;
  for(unsigned int i = 0; i < DS.size(); i++) {
    DF.push_back( ColumnToSEXP(DS.at(i)), DS.colName(i) );
  }
  return DF;
}
