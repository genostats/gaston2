#ifndef SNPWRITEDISK
#define SNPWRITEDISK

#include "SNPvector.h"
#include <vector>
#include <memory> // for shared_ptr
#include "mio.hpp" // for mmap sink
#include <iostream>
#include "mio_utils.h" // for resizing

// TODO : shared ptr ??? on disk only ??? 
// TODO : what do I need ? a boolean if writing went wrong ? 
// TODO : why a class and why inherited ? 
// TODO : combine with the resizing ??????
// TODO : THINK WHY A CLASS
// TODO : PUT FORULA IN A VARIABLE PLEASE

/**
 * @brief Another class derived from SNPvector abstract class (maybe ?)
 * A class to write in a file opened with mio in sink mode
 * @emoji smile
 */
class SNPWriteDisk : public SNPvector {

  public:

  SNPWriteDisk(size_t nbInds, std::shared_ptr<mio::mmap_sink> file_ref, size_t SNP_index, std::vector<uint8_t> to_write, std::string filename) : 
    file_ref_(file_ref), data_( reinterpret_cast<uint8_t *>(file_ref->data() + 3 + (nbInds/4 + ((nbInds%4 == 0u)?0:1)) * SNP_index )) {
    if (!((*file_ref_).is_open())) throw std::runtime_error("File is not registered as properly mapped !\n");
    if ((((*file_ref_).size() - 3) / (nbInds/4 + ((nbInds%4 == 0u)?0:1))) < SNP_index) {
      std::cout << "Index is out of bound, resizing the file to add the SNP at the end...\n";
      if (resizing_file(filename, nbInds/4 + ((nbInds%4 == 0u)?0:1)) < 0) throw std::runtime_error("Resizing of the file failed.\n");
    }
    if (to_write.size() == (nbInds/4 + ((nbInds%4 == 0u)?0:1))) {
      // put size in variable ? 
      size_t starting_index = 3 + (nbInds/4 + ((nbInds%4 == 0u)?0:1)) * SNP_index;
      // to debug :
      std::cout << "File is okay for writing, starting now\n";
      for (size_t i = 0; i < to_write.size(); i++) {
        // cannot use data_
        (*file_ref_).data()[starting_index + i] = to_write[i];
      }
    } else { throw std::runtime_error("The SNP is not of the right size!\n"); }
  }

    // TODO : rm, not necessary
  ~SNPWriteDisk() {
  // std::cout << "Destroying a SNP, here's the count of file_ref_ : " << file_ref_.use_count() << "\n";
  }

  uint8_t * data() {
    return data_;
  }

  /**
   * @brief gives number of indiduals or samples in the SNPs (or else considered an abstract class)
   * 
   * @return size_t 
   */
  size_t nbInds() {
    return nbInds_;
  }
  
  private:
  /** @brief a pointer to the first bit written in the file */
  uint8_t *data_;

   /** @brief to help parse SNP (kept to verify sizof(SNP))*/
  size_t nbInds_;

  size_t SNP_index;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::mmap_sink> file_ref_;
};

#endif // SNPWRITEDISK