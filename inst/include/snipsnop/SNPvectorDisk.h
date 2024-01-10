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

  SNPVectorDisk(size_t nbInds, std::shared_ptr<mio::mmap_source> file_ref, size_t SNP_index) : 
      // data_ is a ptr IN the file, it depends on how many SNPs were read. The 3 offset is to account to the magic nbr in bed file
      data_(file_ref->data() + 3 + (nbInds/4 + ((nbInds%4 == 0u)?0:1)) * SNP_index ), nbInds_(nbInds), file_ref_(file_ref) {} 

  ~SNPVectorDisk() {
  // std::cout << "Destroying a SNP, here's the count of file_ref_ : " << file_ref_.use_count() << "\n";
  }

  size_t nbInds() {
    return nbInds_;
  }

  uint8_t * data() {
    // TODO : think bcos we lost compatibility in reading between SNPvectorMemory
    // TODO : think if not better to have data directly be uint8_t and casted in c°
    // TODO : think of having a way to never go beyond one SNP ( no call more than nbINds/4 + ((nbInds%4 == 0u)?0:1))
    return reinterpret_cast<uint8_t *>(const_cast<char*>(data_));
  }
  
  private:
  /** @brief a pointer to the first bit of SNP in the file */
  const char *data_;
   /** @brief to help parse SNP*/
  size_t nbInds_;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::mmap_source> file_ref_;
};
#endif // SNPMMATRIX