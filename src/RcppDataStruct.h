#include <Rcpp.h>
#include "Column.h"
#include "DataStruct.h"
#include "getDataStructPtr.h"
#include "datatype.h"

#ifndef _rcppdatastruct_
#define _rcppdatastruct_

Rcpp::DataFrame DataStructToDataFrame(const DataStruct & DS);
Rcpp::List DataStructToList(const DataStruct & DS);

inline SEXP ColumnToSEXP(const Column& col) {
  datatype dt = col.type();
  if (dt == INT)
    return Rcpp::wrap(*col.get<int>());
  else if (dt == DOUBLE)
    return Rcpp::wrap(*col.get<double>());
  else if (dt == FLOAT)
    return Rcpp::wrap(*col.get<float>());
  else if (dt == STRING)
    return Rcpp::wrap(*col.get<std::string>());
  else if (dt == BOOL)
    return Rcpp::wrap(*col.get<bool>());
  else
    Rcpp::stop("Unknown column type");
}

#endif
