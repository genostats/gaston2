#include "SNPmatrix.h"
#include "SNPvector.h"
#include "chrType.h"
#include <stdexcept>

// NOTE
// il faut qu'en amont tous les SNPs aient été standardisés de la bonne 
// façon, typiquement par 
// setScaledMode( STANDARDIZED_P, std::sqrt(nbSNPs) );
// pour la GRM classique
//
// NOTE 
// a priori on laissera scalar_t = float même pour remplir une matrice de double
// car ça doit aller plus vite et c'est suffisant
template<typename SNPvectorClass, typename matrixType, typename scalar_t = float>
void GRM(SNPmatrix<SNPvectorClass> & M, matrixType & R) {
  unsigned int nbSNPs = M.nbSNPs();
  unsigned int nbInds = M.nbInds();
  // on vérifie les dimensions de R d'abord
  if(R.nrow() != nbInds || R.ncol() != nbInds)
    throw std::runtime_error("In grm, bad dimensions for R");

  // les résultats bruts [matrice triangulaire, stockée dans un vecteur]
  std::vector<scalar_t> V( (nbInds*(nbInds + 1))/2 );
  for(size_t i = 0; i < nbSNPs; i++) {  
    auto SNP = M.getSNP(i);
    if(SNP->getChrType() == chrType::AUTOSOME) // seulement pour les autosomes
      SNP->template tcrossprod<scalar_t>(V);
  }
 
  // on copie dans R et on symétrise mais ça pourrait être intéressant de
  // renvoyer une matrice symétrique du package Matrix...
  size_t k = 0;
  for(size_t i = 0; i < nbInds; i++) {
    for(size_t j = 0; j <= i; j++) {
      R(j,i) = (double) V[k++];
    }
  }

  // symetriser
  k = 0;
  for(size_t i = 0; i < nbInds; i++) {
    for(size_t j = 0; j <= i; j++) {
      R(i,j) = (double) V[k++]; // ou R(j,i)
    }
  }
}


