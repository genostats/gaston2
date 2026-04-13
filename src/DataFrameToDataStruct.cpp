#include "DataStruct.h"
#include <iostream>
#include <Rcpp.h>


DataStruct DataFrameToDataStruct(Rcpp::DataFrame DF) {
  DataStruct DS;
  SEXP sxpnames = DF.names();
  std::vector<std::string> names;

  if(TYPEOF(sxpnames) == STRSXP) {
    Rcpp::CharacterVector nv(sxpnames);
    for(auto n : nv) {
      names.push_back(Rcpp::as<std::string>(n));
    }
  } else {
    Rcpp::stop("Unnamed data frame ?!");
  }

  for(size_t i = 0; i < DF.size(); i++) {
    SEXP sxp = DF[i];
    switch(TYPEOF(sxp)) {
      case INTSXP: 
        {
          Rcpp::IntegerVector sexv(sxp);
          std::vector<int> v;
          for(int x: sexv) v.push_back(x);
          DS.push_back(v, names[i]);
          break;
        }
      case REALSXP: 
        {
          Rcpp::NumericVector sexv(sxp);
          std::vector<double> v;
          for(double x: sexv) v.push_back(x);
          DS.push_back(v, names[i]);
          break;
        }
      case STRSXP: 
        {
          Rcpp::CharacterVector sexv(sxp);
          std::vector<std::string> v;
          for(auto x: sexv) v.push_back(Rcpp::as<std::string>(x));
          DS.push_back(v, names[i]);
          break;
        }
      default:
        Rcpp::stop("Unsupported column type");
    }
  }

  return DS;
}

