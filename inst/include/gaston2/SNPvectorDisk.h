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


  // constructeur concatenant 2 vecteurs 
  SNPvectorDisk(const std::shared_ptr<SNPvector> first, const std::shared_ptr<SNPvector> second, std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> newfile, size_t SNP_index ) : SNPvector(first->nbInds() + second->nbInds()), file_ref_(newfile) {

    // Initializing properly with the right size and zeros everywhere
    data_ = (uint8_t*)(newfile->data() + 3 + (nbInds_/4 + ((nbInds_%4 == 0u)?0:1)) * SNP_index);
    std::fill(data_, data_ + (nbInds_ / 4 + ((nbInds_ % 4 == 0u) ? 0 : 1)), 0);

    // COPYING THE FIRST SNP
    std::copy(first->data(), first->data() + first->nbChars(), data_);

    // IF the first snp ends right with a byte, I just append it with the second one
    int BitsOnlastByte_first = first->nbInds()%4 * 2;
    if (BitsOnlastByte_first == 0u)  {
      // this should copy everything
      std::copy(second->data(), second->data() + second->nbChars(), data_ + first->nbChars());
    } else {  // FRANKENSTEINING THE SNP :  
      // I cannot simply do that simply because in case of padding I will have holes 
      // so I need to shift everything if that is the case

      uint8_t first_byte = 0;  
      uint8_t second_byte = 0;
      uint8_t carry_over = data_[first->nbChars() - 1];

      int useless_bits_first = 8 - BitsOnlastByte_first;
      int useless_bits_scd = 8 - useless_bits_first;

      // still need to clean up the last byte of first snp just in case
      // to ensure padding is with zeros, and no garbage bits are left for the |= in the loop
      carry_over <<= useless_bits_first;

      size_t new_i = (first->nbChars() - 1);
      // parcourir tout le deuxième à partir du premier 
      for (size_t i = 0; i < second->nbChars(); i++) { // TODO : check que le i soit correct
        // je setup les bytes que je vais manipuler
        first_byte = carry_over;
        second_byte = *(second->data() + i);
        carry_over = second_byte;

        // putting the parts to be merged at the lowest bits possible
        first_byte >>= useless_bits_first;
        data_[new_i] = first_byte;

        // preparing the second byte to be merged on the lowest bits 
        // not conflicting with first byte
        second_byte <<= useless_bits_scd;

        // then mask the first and the second together to write
        data_[new_i] |= second_byte;
        new_i++;
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
