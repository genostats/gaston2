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
                                                             data_(source->data(), source->data() + source->nbChars()) {}

  // constructeur par copie avec une sélection des individus,
  // par un vecteur d'index 
  template <typename intVec>
  SNPvectorMemory(const std::shared_ptr<SNPvector> source, intVec &keep) : SNPvector(keep.size()) {
    const uint8_t *refdata = source->data();
    // check the keep (no idx to far or impossible)
    data_.assign((nbInds_ / 4 + ((nbInds_ % 4 == 0u) ? 0 : 1)), 0); // creating with nbChars of 0

    size_t new_byte = 0;
    size_t new_2bits = 0;

    for (size_t i = 0; i < nbInds_; ++i) {// nbInds_ is already equal to keep.size() here
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
    chr_type_ = source->getChrType();
  }


  // constructeur concatenant 2 vecteurs 
  SNPvectorMemory(const std::shared_ptr<SNPvector> first, const std::shared_ptr<SNPvector> second) : SNPvector(first->nbInds() + second->nbInds()) {

    // Initializing properly with the right size and zeros everywhere
    data_.assign((nbInds_ / 4 + ((nbInds_ % 4 == 0u) ? 0 : 1)), 0); // creating with nbChars of 0

    // COPYING THE FIRST SNP
    //data_(first->data(), first->data() + first->nbChars());
    std::copy(first->data(), first->data() + first->nbChars(), data_.data());

    // IF the first snp ends right with a byte, I just append it with the second one
    int BitsOnlastByte_first = first->nbInds()%4 * 2;
    if (BitsOnlastByte_first == 0u)  {
      // this should copy everything
      std::copy(second->data(), second->data() + second->nbChars(), data_.data() + first->nbChars());
    } else {  // FRANKENSTEINING THE SNP :  
      // I cannot simply do that simply because in case of padding I will have holes 
      // so I need to shift everything if that is the case

      uint8_t first_byte = 0;  
      uint8_t second_byte = 0;
      uint8_t carry_over = data_[first->nbChars() - 1];

      int useless_bits_first = 8 - BitsOnlastByte_first;
      int useless_bits_scd = 8 - useless_bits_first;

      // still need to clean up the last byte of first snp just in case
      // to ensure padding is with zeros, and no garbage bits are left for the |= in the loop
      carry_over <<= useless_bits_first;

      size_t new_i = (first->nbChars() - 1);
      // parcourir tout le deuxième à partir du premier 
      for (size_t i = 0; i < second->nbChars(); i++) { // TODO : check que le i soit correct
        // je setup les bytes que je vais manipuler
        first_byte = carry_over;
        second_byte = *(second->data() + i);
        carry_over = second_byte;

        // putting the parts to be merged at the lowest bits possible
        first_byte >>= useless_bits_first;
        data_[new_i] = first_byte;

        // preparing the second byte to be merged on the lowest bits 
        // not conflicting with first byte
        second_byte <<= useless_bits_scd;

        // then mask the first and the second together to write
        data_[new_i] |= second_byte;
        new_i++;
      }
    }
  }

  uint8_t * data() { return &data_[0]; }

  const uint8_t * data() const { return &data_[0]; }
};

#endif
