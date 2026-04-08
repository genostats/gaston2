#ifndef SNPDOSAGEDISK
#define SNPDOSAGEDISK

#include "SNPdosage.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#include <memory> // for shared_ptr
#include "mio.hpp" // for mmap sink
#include <iostream>

/**
 * @brief Another class derived from SNPdosage abstract class 
 * 
 * Different from SNPdosageMemory because it keeps a shared_ptr
 * to either a mio::mmap_source object or a mio::mmap_sink object
 * depending on the template 'accessMode' parameter (either mio::access_mode::read
 * or mio::access_mode::write), managing the memory mapped file.
 * @emoji smile
 */

template<mio::access_mode accessMode>
class SNPdosageDisk : public SNPdosage {

  private:
  float * data_;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::basic_mmap<accessMode, char>> file_ref_;

  public:
  // /!\ WILL NOT BE ACCOUNTING FOR MAGIC BYTE OFFSET COS CURRENTLY THERE IS NONE 
  //  

  // on donne à ce constructeur un shared ptr vers fichier ouvert par mio + le nb d'individus, et le SNP index à pointer
  SNPdosageDisk(size_t nbInds, std::shared_ptr<mio::basic_mmap<accessMode, char>> file_ref, size_t SNP_index) : 
    SNPdosage(nbInds),
    // no need for sizeof float because the cast is global !
    data_( ((float *) file_ref->data()) + nbInds * SNP_index), 
    file_ref_(file_ref) {}
 

  // TODO !! : check what happens with empty SNP ?
  // constructeur par copie d'un SNPdosage quelconque
  // il va échouer si on n'a pas accessMode == mio::access_mode::write
  // le but est d'écrire dans un fichier .bed *qu'on a créé nous-même* et qui est encore vide (à part les 3 magic bytes)
  // (ou pas forcément vide, ça va la modifier en place)
  // toujours créé en mode RAW_VALUES
  SNPdosageDisk(const std::shared_ptr<SNPdosage> source, 
                std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> file_ref, 
                size_t SNP_index) : 
      SNPdosage(source->nbInds()), 
      data_( ((float *) file_ref->data()) + nbInds_ * SNP_index), 
      file_ref_(file_ref) {
    // on copie les données de source dans le fichier, à la bonne place qui est pointée par data_a
    const float * sourceData = source->data();
    for(size_t i = 0; i < nbInds_; i++) {
      data_[i] = sourceData[i];
    }
  }

  // constructeur par sélection des individus spécifiés : 
  // to_keep contient les indices des individus à conserver dans le SNP
  template <typename intVec>
  SNPdosageDisk(const std::shared_ptr<SNPdosage> source, \
    std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> newfile, size_t SNP_index, intVec &keep) : 
    SNPdosage(keep.size()), file_ref_(newfile) {
    
    data_ = ((float *) newfile->data()) + nbInds_* SNP_index ; // bcos mio sends back a char *
    const float * refdata = source->data();

    for(size_t i = 0; i < nbInds_; i++) {
      int ind_idx = keep[i];
      data_[i] = refdata[ind_idx];
    }
    chr_type_ = source->getChrType();
  }


  // constructeur concatenant 2 vecteurs, et les écrivant dans le nouveau fichier 
  SNPdosageDisk(const std::shared_ptr<SNPdosage> first, const std::shared_ptr<SNPdosage> second, std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> newfile, size_t SNP_index ) 
  : SNPdosage(first->nbInds() + second->nbInds()), file_ref_(newfile) {

    // pointing to the right place for this SNP in the file
    data_ = ((float*) newfile->data()) + nbInds_ * SNP_index;
    const float * firstData = first->data();
    const float * secondData = second->data();

    size_t i = 0;// will use i to write through both snps
    for (; i < first->nbInds(); i++) {
      data_[i] = firstData[i];
    }

    for (size_t j = 0; j < second->nbInds(); j++) {
      data_[i++] = secondData[j];
    }
  }

  ~SNPdosageDisk() {}

  const float * data() const {
    return &data_[0];
  }

};
#endif
