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

template<typename SNPvectorClass, typename matrixType>
void GRM(SNPmatrix<SNPvectorClass> & M, matrixType & R) {
  using scalar_t = typename matrixType::value_type;

  unsigned int nbSNPs = M.nbSNPs();
  unsigned int nbInds = M.nbInds();
  // on vérifie les dimensions de R d'abord
  if(R.nrow() != nbInds || R.ncol() != nbInds)
    throw std::runtime_error("In grm, bad dimensions for R");
  
  // writing in the final matrix BUT as a flat vector
  for(size_t i = 0; i < nbSNPs; i++) {  
    auto SNP = M.getSNP(i);
    if(SNP->getChrType() == chrType::AUTOSOME) // seulement pour les autosomes
      SNP->template tcrossprod<scalar_t>(R);
  }
  
  size_t k = ((nbInds*(nbInds + 1))/2) - 1; //in base 0
  // no need to go through last column...
  for(long j = nbInds -1; j > 0; j--) {
    for(long i = j; i >= 0; i--) {
      R(i,j) = R[k--]; // remplir matrice symétrique à partir de la fin 
      //symétriser icic possiblement : R(j, i) = R(i, j)
    }
  }
  
  // symétriser
  k = 0;
  for(size_t i = 0; i < nbInds; i++) {
    for(size_t j = i + 1; j < nbInds; j++) { // parcours le triangle supérieur
      R(j,i) = R(i,j);
      }
    }
}


