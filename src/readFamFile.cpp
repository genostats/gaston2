#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
void readFamFile(Rcpp::XPtr<SNPmatrix> pM, std::string famFile) {
  pM->readFamFile(famFile);
}

