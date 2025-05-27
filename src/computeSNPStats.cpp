#include "SNPmatrix.h"
#include <Rcpp.h>

// [[Rcpp::export]]
void computeSNPStats(Rcpp::XPtr<SNPmatrix> pM) {
  pM->computeSNPStats();
}
