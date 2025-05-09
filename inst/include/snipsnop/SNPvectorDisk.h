#ifndef SNPVECTORDISK
#define SNPVECTORDISK

#include "SNPvector.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#include <memory> // for shared_ptr
#include "mio.hpp" // for mmap sink
#include <iostream>

/**
 * @brief Another class derived from SNPvector abstract class 
 * 
 * Different from SNPvectorMemory because it keeps a shared_ptr
 * to a mio::mmap_source object managing the memory mapped file.
 * Should the data_ be also a pointer to the data in the file ?
 * @emoji smile
 */
class SNPVectorDisk : public SNPvector {

  public:
  
  SNPVectorDisk(size_t nbInds, std::shared_ptr<mio::mmap_source> file_ref, size_t SNP_index, Mode mode = PLINK) : 
    data_((const uint8_t *) (file_ref->data() + 3 /* offset from the 3 first magic bytes) */ + (nbInds/4 + ((nbInds%4 == 0u)?0:1)) * SNP_index) ), 
    nbInds_(nbInds), file_ref_(file_ref), mode_(mode) {}

  ~SNPVectorDisk() {
    //std::cout << "Destroying a SNP, here's the count of file_ref_ : " << file_ref_.use_count() << "\n";
  }

  size_t nbInds() const { return nbInds_; }

  // pointer to the first char
  // uint8_t * data() {
  //   return &data_[0];
  // }

  const uint8_t * data() const {
    return &data_[0];
  }

  void setMode(Mode mode) {
    mode_ = mode;
  }
  
  void setMode(Mode mode, double personalized[4]) { 
    mode_ = mode;
    for (int i = 0; i< 4; i++)
    currentMode_[4][i] = personalized[i];
  }

  // returns the array used to translate datas
  // TODO : see if Mode enum more usefulS
  const double * mode() const { return currentMode_[mode_]; }
  
  const double mode(unsigned int n) const { return currentMode_[mode_][n]; }
  
  private:
  /** @brief a vector containing the bits composing the SNP */
  //std::vector<uint8_t> data_;
  /** @brief a ptr in the mio file, pointing to the bits composing the SNP
   * It's constness is imposed by mio, that openned a read only file */
  const uint8_t *data_;
   /** @brief to help parse SNP*/
  const size_t nbInds_;
  /** @brief an enum keeping track on how to read datas */
  enum Mode mode_;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::mmap_source> file_ref_;
};
#endif // SNPMMATRIX
