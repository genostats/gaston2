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

  // By default, mode is Numeric (0), TODO : add int checking
  SNPVectorDisk(size_t nbInds, std::shared_ptr<mio::mmap_source> file_ref, int modeInt = 0) : data_(nbInds/4 + ((nbInds%4 == 0u)?0:1)), nbInds_(nbInds), file_ref_(file_ref), mode_((modeInt > 3)? static_cast<SNPvector::Mode>(0) : static_cast<SNPvector::Mode>(modeInt)) {
    if (modeInt > 3 || modeInt < 0) throw std::runtime_error("Wrong mode chosen for reading SNP");
  } 

  ~SNPVectorDisk() {
    //std::cout << "Destroying a SNP, here's the count of file_ref_ : " << file_ref_.use_count() << "\n";
  }

  size_t nbInds() {
    return nbInds_;
  }

  // pointer to the first char
  uint8_t * data() {
    return &data_[0];
  }

  Mode mode() {
    return mode_;
  }
  
  private:
  /** @brief a vector containing the bits composing the SNP */
  std::vector<uint8_t> data_;
  //uint8_t *data_;
   /** @brief to help parse SNP*/
  size_t nbInds_;
  /** @brief an enum helping converting datas if necessary */
  enum Mode mode_;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::mmap_source> file_ref_;
};
#endif // SNPMMATRIX