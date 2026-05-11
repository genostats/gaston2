#include <Rcpp.h>

#include <iostream>

#include "Column.h"
#include "DataStruct.h"
#include "datatype.h"

// [[Rcpp::export]]
void addcolDataStruct(Rcpp::XPtr<DataStruct> pDS, std::string colname, SEXP values) {
  if (!pDS) throw std::runtime_error("data.struct has a broken external pointer !");

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
    default: { // so for single	SINGLESXP and list VECSXP???
      throw std::runtime_error("Added column can only be of type INT, DOUBLE/FLOAT, STRONG or BOOL");
    }
  }
}