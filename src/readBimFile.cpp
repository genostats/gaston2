#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
void readBimFile(Rcpp::XPtr<SNPmatrix> pM, std::string bimFile) {
  pM->readBimFile(bimFile);
}

