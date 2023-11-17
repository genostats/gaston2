#include "SNPvector.h"
#include <cstddef>
#include <cstdint>
#include <vector>

#ifndef _snpvectormemory_
#define _snpvectormemory_
class SNPvectorMemory : public SNPvector {
  private:
  std::vector<uint8_t> data_;
  size_t nbInds_;
  public:
  // un constructeur, taille nbInds
  SNPvectorMemory(size_t nbInds) : data_(nbInds/4 + ((nbInds%4 == 0u)?0:1)), nbInds_(nbInds) {} 


  // un constructeur à partir d'un vecteur de char 
  // !! ce constructeur fait un move et donc "vide" son argument
  // !! est-ce qu'on en a besoin ?! -> commenté pour le moment
  /*
  SNPvectorMemory(std::vector<uint8_t> & data) : data_(std::move(data)) {}
  */

  // nb of individuals
  size_t nbInds() {
    return nbInds_;
  }
  // pointer to the first char
  uint8_t * data() {
    return &data_[0];
  } 
};

#endif
