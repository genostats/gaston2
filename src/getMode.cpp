#include "SNPmatrix.h"
#include "mode.h"
#include <Rcpp.h>


// [[Rcpp::export]]
std::string getMode(Rcpp::XPtr<SNPmatrix<>> pM) {
  return modeToString(pM->mode());
}

