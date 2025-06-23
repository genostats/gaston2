#include "SNPmatrix.h"
#include "GRM.h"
#include "chrType.h"
#include <iostream>
#include <stdexcept>
#include <Rcpp.h>
#include "debug.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix grm(Rcpp::XPtr<SNPmatrix<>> pM) {
  unsigned int nbSNPs = pM->nbSNPs();
  unsigned int nbInds = pM->nbInds();
  Rcpp::NumericMatrix R(nbInds, nbInds);

  // TODO repenser à ces préliminaires,
  // est-ce que c'est souhaitable d'en faire un membre de SNPmatrix ?

  // -------- préparation appel GRM ----------
  //
  // il faut que tous les SNPs soient 
  // convenablement standardisés
  
  // d'abord on calcule les stats
  pM->computeSNPStats();

  // on compte les SNPs autosomaux
  int nbAutosomalSNPs = 0;
  for(size_t i = 0; i < nbSNPs; i++) {  
    auto SNP = pM->getSNP(i);
    if(SNP->getChrType() == chrType::AUTOSOME) nbAutosomalSNPs++;
  }

  if(nbAutosomalSNPs == 0) Rcpp::stop("No autosomal SNPs");

  // on les met à l'échelle voulue
  for(size_t i = 0; i < nbSNPs; i++) {  
    auto SNP = pM->getSNP(i);
    SNP->setScaledMode( STANDARDIZED_P, std::sqrt(nbAutosomalSNPs) );
  }
  // ----------------------------------------------

  // calcul effectif de la GRM
  GRM(*pM, R);  

  return R;
}


