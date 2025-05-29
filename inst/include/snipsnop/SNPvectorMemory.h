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
class SNPvectorMemory : public SNPvector
{
public:
  // un constructeur, taille nbInds
  SNPvectorMemory(size_t nbInds) : data_(nbInds / 4 + ((nbInds % 4 == 0u) ? 0 : 1)), nbInds_(nbInds) {}

  // constructeur par copie
  SNPvectorMemory(const std::shared_ptr<SNPvector> source) : data_(source->data(), source->data() + source->nbInds()),
                                                             nbInds_(source->nbInds()) {}

  // constructeur par copie avec une sélection des individus,
  // par un vecteur d'index 
  template <typename intVec>
  SNPvectorMemory(const std::shared_ptr<SNPvector> source, intVec &keep) : nbInds_(keep.size()) {
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


  // nb of individuals
  size_t nbInds() const { return nbInds_; }

  /**
   * @brief A function pointing to the first character of the SNP,
   * It's the only class from SNPVector to have a non_const data implemented.
   * It is for writing data to the freshly new snp.
   *
   * @return uint8_t*
   */
  uint8_t *data() { return &data_[0]; }

  const uint8_t *data() const { return &data_[0]; }

  // void setMode(Mode mode) { mode_ = mode; }
  // void setMode(Mode mode, double personalized[4]) { 
  //  mode_ = mode;
  //  for (int i = 0; i< 4; i++)
  //   currentMode_[4][i] = personalized[i];
  // }

  // returns the array used to translate datas 
  // const double *mode() const { return currentMode_[mode_]; } // considered as const cos double *mode() const

  // returns the translation of n through mode
  // const double mode(unsigned int n) const { return currentMode_[mode_][n]; }

  private:
  std::vector<uint8_t> data_;
  const size_t nbInds_;

};

#endif
