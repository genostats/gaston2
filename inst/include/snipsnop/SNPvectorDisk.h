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

  public:
 
  // on donne à ce constructeur un shared ptr vers fichier ouvert par mio + le nb d'individus, et le SNP index à pointer
  SNPvectorDisk(size_t nbInds, std::shared_ptr<mio::basic_mmap<accessMode, char>> file_ref, size_t SNP_index, Mode mode = PLINK) : 
    data_((uint8_t *) (file_ref->data() + 3 /* offset from the 3 first magic bytes) */ + (nbInds/4 + ((nbInds%4 == 0u)?0:1)) * SNP_index) ), 
    nbInds_(nbInds), file_ref_(file_ref), mode_(mode) {}

  // TODO !! : check what happens with empty SNP ? 

  // constructeur par copie d'un SNPVector quelconque
  // il va échouer si on n'a pas accessMode == mio::access_mode::write
  // le but est d'écrire dans un fichier .bed *qu'on a créé nous-même* et qui est encore vide (à part les 3 magic bytes)
  // (ou pas forcément vide, ça va la modifier en place)
  // toujours créé en mode PLINK
  SNPvectorDisk(const std::shared_ptr<SNPvector> source, std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> file_ref, size_t SNP_index) : 
      data_((uint8_t *) (file_ref->data() + 3 + (source->nbInds()/4 + ((source->nbInds()%4 == 0u)?0:1)) * SNP_index)), nbInds_(source->nbInds()), file_ref_(file_ref), mode_(PLINK) {
    // on copie les données de source dans le fichier, à la bonne place qui est pointée par data_a
    size_t nbChars = source->nbChars();
    const uint8_t * sourceData = source->data();
    for(size_t i = 0; i < nbChars; i++) {
      data_[i] = sourceData[i];
    }
  }

  // constructeur par sélection des individus spécifiés : 
  // to_keep contient les indexs des individus à conserver dans le SNP
  SNPvectorDisk(const std::shared_ptr<SNPvector> source, std::shared_ptr<mio::basic_mmap<mio::access_mode::write, char>> newfile, size_t SNP_index, const std::vector<size_t> &keep) : 
    nbInds_(keep.size()), file_ref_(newfile), mode_(PLINK) {
    
    data_ = (uint8_t *) newfile->data() + (3 + (nbInds_/4 + ((nbInds_%4 == 0u)?0:1)) * SNP_index); // bcos mio sends back a char *
    const uint8_t * sourceData = source->data();
    uint8_t newdata = 0;

    for(size_t i = 0; i < nbInds_; i++) {

      int ind_idx = keep[i];
      size_t currentChar = ind_idx / 4;    // index du byte
      int ind_gen = read_ind(sourceData[currentChar], ind_idx);

      size_t new_byte = (i / 4);
      size_t new_2bits = (i % 4) * 2;

      ind_gen <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte
      newdata |= ind_gen;

      //size_t snp_offset = SNP_index * ((nbInds_ + 3) / 4); // ensure rounding up
      /* If i want to write byte by byte instead of 2bits by 2bits, this is the way*/
      // once every 4 inds BYTE BY BYTE
      if (i % 4 == 3 || i == nbInds_ - 1)
      {
        data_[new_byte] = newdata;
        newdata = 0;
      }
    }
  }

  ~SNPvectorDisk() {
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
  uint8_t *data_;
   /** @brief to help parse SNP*/
  const size_t nbInds_;
  /** @brief an enum keeping track on how to read datas */
  enum Mode mode_;
  /** @brief a shared_ptr to the object handling the file, 
  when the last SNP from the same file is deleted, the file is unmapped and closed */
  std::shared_ptr<mio::basic_mmap<accessMode, char>> file_ref_;
};
#endif // SNPMMATRIX
