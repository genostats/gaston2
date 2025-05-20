#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerMatrix SNPMatrixToIntegerMatrix(Rcpp::XPtr<SNPmatrix> pM) {
  unsigned int nbSNPs = pM->nbSNPs();
  unsigned int nbInds = pM->nbInds();
  Rcpp::IntegerMatrix R(nbInds, nbSNPs);
  for(unsigned int j = 0; j < nbSNPs; j++) {
    auto SNP = pM->getSNP(j);
    SNP->setMode(PLINK);
    unsigned int i = 0;
    for(int x : *SNP) {
      R(i++, j) = (x < 3)?x:NA_INTEGER;
    }
  }
  return R;
}


