#include <Rcpp.h>

#include "SNPmatrix.h"
#include "GRM.h"
#include "chrType.h" // pour le check sur les autosomes
#include <stdexcept>

#include "houba/MMatrix.h"

// [[Rcpp::export]]
SEXP grm_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, std::string path) {
  unsigned int nbSNPs = pM->nbSNPs();
  unsigned int nbInds = pM->nbInds();
  
  // -------- préparation objet houba ----------

  Rcpp::S4 GRM_MMatrix("mmatrix");
  GRM_MMatrix.slot("file") = path;
  GRM_MMatrix.slot("dim") = Rcpp::IntegerVector::create(nbInds, nbInds);
  GRM_MMatrix.slot("readonly") =  false;

  // cast en double au moment du calcul dans GRM.h
  Rcpp::XPtr<houba::MMatrix<double>>  GRM_ptr(new houba::MMatrix<double>(path, nbInds, nbInds));
  // make sure that the file is zeroed (if it exists before hand, its bytes are still there)
  size_t matrix_size = nbInds * nbInds * sizeof(double);
  std::memset(GRM_ptr->data(), 0, matrix_size);

  GRM_MMatrix.slot("ptr") = GRM_ptr;
  GRM_MMatrix.slot("datatype") = "double";

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
    if(SNP->getChrType() == chrType::AUTOSOME && !SNP->monomorphe()) nbAutosomalSNPs++;
  }

  if(nbAutosomalSNPs == 0) Rcpp::stop("No autosomal (non monomorphic) SNPs");

  // on les met à l'échelle voulue
  for(size_t i = 0; i < nbSNPs; i++) {  
    auto SNP = pM->getSNP(i);
    SNP->setScaledMode( STANDARDIZED_P, std::sqrt(nbAutosomalSNPs) );
  }
  // ----------------------------------------------

  // calcul effectif de la GRM

  GRM(*pM, *GRM_ptr);  

  //std::cout << "Message from houba : " << GRM_ptr->getVerbosout() << "\n"; // TODO REMOVE, just to debug

  return GRM_MMatrix;
}