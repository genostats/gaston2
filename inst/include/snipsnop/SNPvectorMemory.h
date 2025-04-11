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
  SNPvectorMemory(size_t nbInds, int modeInt = 0) : data_(nbInds / 4 + ((nbInds % 4 == 0u) ? 0 : 1)), nbInds_(nbInds), mode_((modeInt > 3) ? static_cast<Mode>(0) : static_cast<Mode>(modeInt))
  {
    if (modeInt > 3 || modeInt < 0)
      throw std::runtime_error("Wrong mode chosen for reading SNP");
  }

  // un constructeur à partir d'un vecteur de char
  // !! ce constructeur fait un move et donc "vide" son argument
  // !! est-ce qu'on en a besoin ?! -> commenté pour le moment
  /*
  SNPvectorMemory(std::vector<uint8_t> & data) : data_(std::move(data)) {}
  */

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

  void setMode(Mode mode) { mode_ = mode; }

  // returns the array used to translate datas 
  const double *mode() const { return currentMode_[mode_]; } // considered as const cos double *mode() const

  //returns the translation of n through mode
  const double mode(unsigned int n) const { return currentMode_[mode_][n]; }

  private:
  std::vector<uint8_t> data_;
  const size_t nbInds_;
  enum Mode mode_;

};

#endif