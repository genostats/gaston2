#include "DataStruct.h"
#include "SNPmatrix.h"
#include <Rcpp.h>

DataStruct * getDataStructPtr(Rcpp::S4 x) {
  int type = x.slot("type");
  DataStruct * pDS;
  switch(type) {
    case 0: {
      pDS = Rcpp::as<Rcpp::XPtr<DataStruct>>(x.slot("ptr"));
      break;
    }
    case 1: {
      pDS = &Rcpp::as<Rcpp::XPtr<SNPmatrix<>>>(x.slot("ptr"))->getIndStats();
      break;
    }
    case 2: {
      pDS = &Rcpp::as<Rcpp::XPtr<SNPmatrix<>>>(x.slot("ptr"))->getSNPStats();
      break; 
    }
    default:
      throw std::runtime_error("unknown type");
  }
  if(!pDS) throw std::runtime_error("unvalid ptr");
  return pDS;
}

