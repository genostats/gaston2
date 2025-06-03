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

  public:
  // un constructeur, taille nbInds
  SNPvectorMemory(size_t nbInds) : SNPvector(nbInds), data_(nbInds / 4 + ((nbInds % 4 == 0u) ? 0 : 1)) {}

  // constructeur par copie
  SNPvectorMemory(const std::shared_ptr<SNPvector> source) : SNPvector(source->nbInds()), 
                                                             data_(source->data(), source->data() + source->nbInds()) {}

  // constructeur par copie avec une sélection des individus,
  // par un vecteur d'index 
  template <typename intVec>
  SNPvectorMemory(const std::shared_ptr<SNPvector> source, intVec &keep) : SNPvector(keep.size()) {
    const uint8_t *refdata = source->data();
    // check the keep (no idx to far or impossible)
    data_.assign((nbInds_ / 4 + ((nbInds_ % 4 == 0u) ? 0 : 1)), 0); // creating with nbChars of 0

    size_t new_byte = 0;
    size_t new_2bits = 0;

    for (size_t i = 0; i < nbInds_; ++i) {
      int ind_idx = keep[i];
      size_t currentChar = ind_idx / 4;
      int ind_gen = read_ind(refdata[currentChar], ind_idx);
      
      ind_gen <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte

      new_2bits += 2;
      data_[new_byte] |= ind_gen; // le || pour ne pas toucher au bits déjà set
      if (new_2bits == 8){
        new_byte++;
        new_2bits = 0;
      }
    }
  }

  uint8_t * data() { return &data_[0]; }

  const uint8_t * data() const { return &data_[0]; }
};

#endif
