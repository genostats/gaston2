#include "SNPmatrix.h"
#include <Rcpp.h>

// [[Rcpp::export]]
void setChrType(Rcpp::XPtr<SNPmatrix<SNPvector>> pM) {
  pM->setChrType();
}
