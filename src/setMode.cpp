#include "SNPmatrix.h"
#include "mode.h"
#include <Rcpp.h>


// [[Rcpp::export]]
void setMode(Rcpp::XPtr<SNPmatrix<>> pM, std::string mode) {
  pM->computeSNPStats();
  pM->setMode( stringToMode(mode) );
}

