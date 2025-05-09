#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerMatrix SNPMatrixToIntegerMatrix(Rcpp::XPtr<SNPmatrix> pM) {
  unsigned int nbSNPs = pM->nSNP();
  if(nbSNPs == 0) {
    Rcpp::IntegerMatrix R(0, 0);
    return R;
  }
  auto SNP0 = pM->getSNP(0);
  unsigned int nbInds = SNP0->nbInds();
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


