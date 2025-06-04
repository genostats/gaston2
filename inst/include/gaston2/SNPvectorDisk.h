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
 * to either a mio::mmap_source object or a mio::mmap_sink object
 * depending on the template 'accessMode' parameter (either mio::access_mode::read
 * or mio::access_mode::write), managing the memory mapped file.
 * @emoji smile
 */

template<mio::access_mode accessMode>
class SNPvectorDisk : public SNPvector {

  private:
  uint8_t * data_;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::basic_mmap<accessMode, char>> file_ref_;

  public:
  // on donne à ce constructeur un shared ptr vers fichier ouvert par mio + le nb d'individus, et le SNP index à pointer
  SNPvectorDisk(size_t nbInds, std::shared_ptr<mio::basic_mmap<accessMode, char>> file_ref, size_t SNP_index) : 
    SNPvector(nbInds),
    data_((uint8_t *) (file_ref->data() + 3 /* offset from the 3 first magic bytes) */ + (nbInds/4 + ((nbInds%4 == 0u)?0:1)) * SNP_index) ), 
    file_ref_(file_ref) {}
 

  // TODO !! : check what happens with empty SNP ?
  // constructeur par copie d'un SNPVector quelconque
  // il va échouer si on n'a pas accessMode == mio::access_mode::write
  // le but est d'écrire dans un fichier .bed *qu'on a créé nous-même* et qui est encore vide (à part les 3 magic bytes)
  // (ou pas forcément vide, ça va la modifier en place)
  // toujours créé en mode RAW_VALUES
  SNPvectorDisk(const std::shared_ptr<SNPvector> source, 
                std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> file_ref, 
                size_t SNP_index) : 
      SNPvector(source->nbInds()), 
      data_((uint8_t *) (file_ref->data() + 3 + (source->nbInds()/4 + ((source->nbInds()%4 == 0u)?0:1)) * SNP_index)), 
      file_ref_(file_ref) {
    // on copie les données de source dans le fichier, à la bonne place qui est pointée par data_a
    size_t nbChars = source->nbChars();
    const uint8_t * sourceData = source->data();
    for(size_t i = 0; i < nbChars; i++) {
      data_[i] = sourceData[i];
    }
  }

  // constructeur par sélection des individus spécifiés : 
  // to_keep contient les indices des individus à conserver dans le SNP
  template <typename intVec>
  SNPvectorDisk(const std::shared_ptr<SNPvector> source, \
    std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> newfile, size_t SNP_index, intVec &keep) : 
    SNPvector(keep.size()), file_ref_(newfile) {
    
    data_ = (uint8_t *) newfile->data() + (3 + (nbInds_/4 + ((nbInds_%4 == 0u)?0:1)) * SNP_index); // bcos mio sends back a char *
    const uint8_t * sourceData = source->data();
    uint8_t newdata = 0;
    size_t new_byte = 0;
    size_t new_2bits = 0;

    for(size_t i = 0; i < nbInds_; i++) {

      size_t ind_idx = keep[i];
      size_t currentChar = ind_idx / 4;    // index du byte
      uint8_t ind_gen = read_ind(sourceData[currentChar], ind_idx);

      ind_gen <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte
      newdata |= ind_gen;

      new_2bits += 2;
      /* If i want to write byte by byte instead of 2bits by 2bits, this is the way*/
      if (new_2bits == 8 || i == nbInds_ - 1)
      {
        data_[new_byte++] = newdata;
        newdata = 0;
        new_2bits = 0;
      }
    }
  }

  ~SNPvectorDisk() {
    //std::cout << "Destroying a SNP, here's the count of file_ref_ : " << file_ref_.use_count() << "\n";
  }

  const uint8_t * data() const {
    return &data_[0];
  }

};
#endif // SNPMMATRIX
