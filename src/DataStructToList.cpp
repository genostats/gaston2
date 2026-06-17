#include <Rcpp.h>
#include "RcppDataStruct.h"

Rcpp::List DataStructToList(const DataStruct & DS) {
  Rcpp::List L;
  for(unsigned int i = 0; i < DS.size(); i++) {
    L.push_back( ColumnToSEXP(DS.at(i)), DS.colName(i) );
  }
  return L;
}

