#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix SNPMatrixToNumericMatrix(Rcpp::XPtr<SNPmatrix<>> pM) {
  unsigned int nbSNPs = pM->nbSNPs();
  unsigned int nbInds = pM->nbInds();
  Rcpp::NumericMatrix R(nbInds, nbSNPs);
  for(unsigned int j = 0; j < nbSNPs; j++) {
    auto SNP = pM->getSNP(j);
    SNP->compute_stats();
    SNP->setMode(STANDARDIZED_P);
    unsigned int i = 0;
    for(double x : *SNP) {
      R(i++, j) = x;
    }
  }
  return R;
}


