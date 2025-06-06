#include "SNPmatrix.h"
#include <Rcpp.h>

// [[Rcpp::export]]
void computeSNPStats(Rcpp::XPtr<SNPmatrix<SNPvector>> pM) {
  pM->computeSNPStats();
}
