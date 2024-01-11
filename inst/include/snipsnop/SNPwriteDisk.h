#ifndef SNPWRITEDISK
#define SNPWRITEDISK



// TODO : alleviate includes pls 
#include "SNPvector.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#include <memory> // for shared_ptr
#include "mio.hpp" // for mmap sink
#include <iostream>

// TODO : shared ptr ??? on disk only ??? 
// TODO : what do I need ? a file to write to, an index, and a data, a boolean if writing went wrong ? 
// TODO : why a class and why inherited ? 
// TODO : write in insertion ? by replacing ? 
// TODO : combine with the resizing ??????
// TODO : THINK WHY A CLASS

// TODO : think, by doing that need to update current SNP 
/**
 * @brief Another class derived from SNPvector abstract class (maybe ?)
 * A class to write in a file opened with mio in sink mode
 * @emoji smile
 */
class SNPWriteDisk : public SNPvector {

  public:

  SNPWriteDisk(size_t nbInds, std::shared_ptr<mio::mmap_sink> file_ref, size_t SNP_index, std::vector<uint8_t> to_write) : 
    file_ref_(file_ref), data_( reinterpret_cast<uint8_t *>(file_ref->data() + 3 + (nbInds/4 + ((nbInds%4 == 0u)?0:1)) * SNP_index )) {
    if (!((*file_ref_).is_open())) std::cerr << "File is not registered as properly mapped !\n";
    if (((*file_ref_).size() - 3) / 8 < SNP_index) std::cerr << "Index seems out of bound\n"; // eventually resizing it + CHECK MODULO
    if (to_write.size() == (nbInds/4 + ((nbInds%4 == 0u)?0:1))) {
      // write, so for boucle on every bit
      // put size in variable ? 
      // to debug :
      std::cout << "File is okay for writing, starting now\n";
      for (int i = 0; i < to_write.size(); i++) {
        data_[i] = to_write[i];
      }
    } else {
      //error, but exception ? std::cerr ? 
    }
    


  } 

    // TODO : rm, not necessary
  ~SNPWriteDisk() {
  // std::cout << "Destroying a SNP, here's the count of file_ref_ : " << file_ref_.use_count() << "\n";
  }

  uint8_t * data() {
    return data_;
  }
  
  private:
  /** @brief a pointer to the first bit where to write in the file */
  uint8_t *data_;

   /** @brief to help parse SNP (kept to verify sizof(SNP))*/
  size_t nbInds_;

  size_t SNP_index;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::mmap_sink> file_ref_;
};


#endif // SNPWRITEDISK