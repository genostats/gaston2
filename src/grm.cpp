#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix grm(Rcpp::XPtr<SNPmatrix> pM) {
  unsigned int nbSNPs = pM->nbSNPs();
  unsigned int nbInds = pM->nbInds();
  std::vector<float> V( (nbInds*(nbInds + 1))/2 );
  pM->computeSNPStats();
  for(size_t i = 0; i < nbSNPs; i++) {  
    auto SNP = pM->getSNP(i);
    // TODO repenser à ceci... peut-etre créer une fonction dans SNPmatrix.h pour
    // appliquer à tous les SNPs ?
    SNP->setScaledMode( STANDARDIZED_P, std::sqrt(nbSNPs) );
    SNP->tcrossprod<float>(V);
  }
 
  // on symmétrise mais ça pourrait être intéressant de
  // renvoyer une matrice symétrique du package Matrix...!
  Rcpp::NumericMatrix R(nbInds, nbInds);
  size_t k = 0;
  for(size_t i = 0; i < nbInds; i++) {
    for(size_t j = 0; j <= i; j++) {
      R(j,i) = (double) V[k++];
    }
  }

  // symmetriser
  k = 0;
  for(size_t i = 0; i < nbInds; i++) {
    for(size_t j = 0; j <= i; j++) {
      R(i,j) = (double) V[k++]; // ou R(j,i)
    }
  }
  return R;
}


