#include "SNPvector.h"
#include <cstddef>
#include <cstdint>
#include <vector>

#ifndef _snpvectormemory_
#define _snpvectormemory_


/**
 * @brief A class derived from SNPvector abstract class 
 * 
 * keeps SNP data in memory and 
 * @emoji cry
 */
class SNPvectorMemory : public SNPvector {
  private:
  std::vector<uint8_t> data_;
  size_t nbInds_;
  const uint8_t* mode_;
  public:
  // un constructeur, taille nbInds
  SNPvectorMemory(size_t nbInds, const uint8_t* mode = Defaultmode ) : data_(nbInds/4 + ((nbInds%4 == 0u)?0:1)), nbInds_(nbInds), mode_(mode) {
    if (!mode) {
        mode_ = Defaultmode;
    }
  } 


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

  /**
   * @brief A function pointing to the first character of the SNP 
   * @example 
   * 01101100
   * ^ here 
   * 
   * @return uint8_t* 
   */
  uint8_t * data() {
    return &data_[0];
  } 

  const uint8_t * mode() {
    return mode_;
  }

  int mode(unsigned int n) {
    return mode_[n];
  }
};

#endif