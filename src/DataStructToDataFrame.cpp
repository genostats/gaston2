#include <Rcpp.h>

#include <iostream>

#include "DataStruct.h"

SEXP ColumnToSEXP(const Column& col) {
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


// When column not all the same size, will downgrade to a Rcpp::List
Rcpp::DataFrame DataStructToDataFrame(const DataStruct& DS) {
  size_t nb_col = DS.size();
  // First need to run a check of if all column have same size
  size_t ref_size = DS.at(0).size();
  unsigned int i = 0;
  while (i < nb_col) {
    if (DS.at(i).size() != ref_size) break;
    i++;
  }

  if (i == nb_col - 1) {
    // all same size, DF possible
    Rcpp::DataFrame DF;
    for (unsigned int i = 0; i < nb_col; i++) {
      DF.push_back(ColumnToSEXP(DS.at(i)), DS.colName(i));
    }
    return DF;
  } else {
    // downgrade to listo remove 
    // TODO :to rm debug
    std::cout << "using OUR list \n";
    Rcpp::List L;
    for (unsigned int i = 0; i < nb_col; i++) {
      L.push_back(ColumnToSEXP(DS.at(i)), DS.colName(i));
    }
    return L;
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame DataStructToDataFrame_(Rcpp::XPtr<DataStruct> pDS) {
  return DataStructToDataFrame(*pDS);
}