#include <Rcpp.h>
#include "RcppDataStruct.h"

Rcpp::List DataStructToList(const DataStruct & DS, Rcpp::IntegerVector I) {
  Rcpp::List L;
  for(int i : I) {
    L.push_back( ColumnToSEXP(DS.at(i)), DS.colName(i) );
  }
  return L;
}

