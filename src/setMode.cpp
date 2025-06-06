#include "SNPmatrix.h"
#include "mode.h"
#include <Rcpp.h>


// [[Rcpp::export]]
void setMode(Rcpp::XPtr<SNPmatrix<>> pM, std::string mode) {
  pM->setMode( stringToMode(mode) );
}

