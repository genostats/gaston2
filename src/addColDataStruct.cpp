#include <Rcpp.h>

#include <iostream>

#include "Column.h"
#include "DataStruct.h"
#include "getDataStructPtr.h"
#include "datatype.h"

// [[Rcpp::export]]
void addColDataStruct(Rcpp::S4 x, std::string colname, SEXP values) {
  DataStruct * pDS = getDataStructPtr(x);

  // must create a Column object
  // then Call pDS->setColumn(const Column & col, std::string name) {}
  // will check if it already exists and potentially replace it ! Else will just add it
  switch (TYPEOF(values)) {
    case INTSXP: { // int
      return pDS->setColumn(Column(Rcpp::as<std::vector<int>>(values)) ,colname);
    }
    case REALSXP: { // double AND float TODO : think of a solution for that
      return pDS->setColumn(Column(Rcpp::as<std::vector<double>>(values)) ,colname);
    }
    case STRSXP: { // string
      return pDS->setColumn(Column(Rcpp::as<std::vector<std::string>>(values)) ,colname);
    }
    case LGLSXP: {
      return pDS->setColumn(Column(Rcpp::as<std::vector<bool>>(values)), colname);
    }
    default: { 
      throw std::runtime_error("Added column can only be of type INT, DOUBLE/FLOAT, STRONG or BOOL");
    }
  }
}
