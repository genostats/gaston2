#ifndef SNPWRITEDISK
#define SNPWRITEDISK

#include <vector>
#include <memory> // for shared_ptr
#include "mio.hpp" // for mmap sink
#include <iostream>
#include "mio_utils.h" // for resizing

// TODO : what do I need ? a boolean if writing went wrong ? 
// TODO : why a class and why inherited ? can keep the shared ptr and so the file opened eg
// TODO: check that nothing went wrong with the constructor, and throwing good way of handling errorsd

/**
 * @brief Another class derived from SNPvector abstract class (maybe ?)
 * A class to write in a file opened with mio in sink mode
 * @emoji smile
 */
class SNPWriteDisk {

  public:

  SNPWriteDisk(size_t nbInds, std::shared_ptr<mio::mmap_sink> file_ref, size_t SNP_index, std::vector<uint8_t> to_write, std::string filename) : 
    file_ref_(file_ref), sizeofSNP_(nbInds/4 + ((nbInds%4 == 0u)?0:1))  {
    if (!((*file_ref_).is_open())) throw std::runtime_error("File is not registered as properly mapped !\n");
    if ((((*file_ref_).size() - 3) / sizeofSNP_) < SNP_index) {
      std::cout << "Index is out of bound, resizing the file to add the SNP at the end...\n";
      if (resizing_file(filename, sizeofSNP_) < 0) throw std::runtime_error("Resizing of the file failed.\n");
    }
    if (to_write.size() == sizeofSNP_) {
      size_t starting_index = 3 + (sizeofSNP_ * SNP_index);
      // to debug :
      std::cout << "File is okay for writing, starting now\n";
      for (size_t i = 0; i < to_write.size(); i++) {
        // cannot use data_
        (*file_ref_).data()[starting_index + i] = to_write[i];
      }
    } else { throw std::runtime_error("The SNP is not of the right size!\n"); }
  }
  
  private:
  size_t sizeofSNP_;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::mmap_sink> file_ref_;
};

#endif // SNPWRITEDISK