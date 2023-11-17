#include <cstdint>
#include <cstddef>
#ifndef _SNPvector_
#define _SNPvector_

class SNPvector {
  public:
  // fonctions purement virtuelles : 
  // récupérer un pointeur vers le début du vecteur d'unsigned chars
  virtual uint8_t * data() = 0;
  // nombre d'individus
  virtual size_t nbInds() = 0;
 
  // size of the vector
  size_t nbChars() {
    size_t n = nbInds();
    return n/4 + ((n%4 == 0u)?0:1);
  }

  // an idiotic function just for testing purposes
  unsigned int sum() {
    unsigned int S = 0;
    uint8_t * d = data();
    size_t nc = nbChars();
    for(size_t i = 0; i < nc; i++) {
      S += (unsigned int) d[i];
    }
    return S;
  }
};

#endif
