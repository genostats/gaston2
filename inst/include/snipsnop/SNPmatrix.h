#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector>
#include <memory>
#include <array>     // for compute stats function
#include <stdexcept> // for out of range exceptions

#include "SNPvectorMemory.h"

#include <omp.h>

#ifndef _snpmatrix_
#define _snpmatrix_

/**
 * @brief A class keeping in a vector pointers to SNPvectors
 * (either SNPvectorDisk or SNPvectorMemory)
 *
 */
class SNPmatrix
{
public:
    /**
   * @brief function to add a shared_ptr into the vector keeping the 
   * SNPs of the matrix. /!\ will check if the SNP added is matching in size
   *
   */
  void push_back(std::shared_ptr<SNPvector> v)
  {
    // if at least one loaded, everySNPs must have same size
    if (SNPs.size() && (SNPs[0]->nbInds() != v->nbInds()))
    {
      std::cerr << "Pb loading SNP" << std::endl;
      throw std::out_of_range("Attempting to load a SNP with a different nb of individuals");
    }
    SNPs.push_back(v);
  }

  int size() { return SNPs.size(); }

  // temporary func to test d°
  void deleteSNP()
  {
    SNPs.pop_back();
  }

  /**
   * @brief A member function calculating N0s, N1s, N2s and NAs (in Plink format) for one individual
   * in the SNPmatrix.
   *
   * @param Indnb the individual to look for in every snps
   * @return std::array<int, 4> with [0] = N0s, [1] = N1s, [2] = N2s and [3] = NAs
   */
  std::array<int, 4> compute_stats(int Indnb)
  {
    int total = SNPs.size();
    if (total == 0)
      throw std::out_of_range("No SNPs loaded into the SNPMatrix !");
    if (Indnb > SNPs[0]->nbInds() - 1)
      throw std::out_of_range("You're trying to access an out-of bound indiviual.");
    std::array<int, 4> stats = {{0, 0, 0, 0}};

#pragma omp parallel for
    for (int i = 0; i < total; i++)
    { // parcourt tous les SNPs
      auto vec = SNPs[i];
      int j = 0;
      for (auto pa = vec->begin(); pa != vec->end(); ++pa)
      { // parcourt le SNP[i]
        if (j == Indnb)
        {
/* is not good practice to put atomic inside the loop,
 considered "serializing your code",
means every thread will synchronize for this
+ not a scalar (array), so false
but I guess it's fine when after a condition (relatively) rarely met ? */
#pragma omp atomic
          stats[(int)*pa]++;
          break;
        }
        j++;
      }
    }
    return stats;
  }

  // WIP : NEED TO TEST THOROUGHLY, FIND A WAY FOR SNPVECTORDISK
  SNPmatrix extract_ind(const std::vector<int> &keep)
  {

    if (SNPs.empty())
      throw std::runtime_error("Original matrix is empty !");

    size_t ref_inds = SNPs[0]->nbInds();
    for (int i : keep)
    {
      if (i < 0 || i >= ref_inds)
        throw std::out_of_range("No individuals with this index : " + std::to_string(i));
    }

    SNPmatrix new_matrix;
    size_t new_nb_inds = keep.size();

      for (const auto &snp : SNPs)
      { // parcours de tous les snps de la matrices d'origine

        size_t i = 0;
        if (!this->onDisk(i++)) {
        // to feed to new matrix
        std::shared_ptr<SNPvectorMemory> new_snp = std::make_shared<SNPvectorMemory>(new_nb_inds);
        const uint8_t *refdata = snp->data();
          
        // Filling with 0, #https://www.cog-genomics.org/plink/1.9/formats#bed
        // If N is not divisible by four, the extra high-order bits in the last byte of each block are always zero.
        std::fill(new_snp->data(), new_snp->data() + new_snp->nbChars(), 0);

        for (size_t i = 0; i < new_nb_inds; ++i)
        { // parcours de keep pour isoler lees individus
          int ind_idx = keep[i];
          // recherche de l'individus
          size_t currentChar = ind_idx / 4;    // index du byte
          uint8_t byte = refdata[currentChar]; // je récup la byte
          size_t current2bits = (ind_idx % 4); // where bits dans le byte
          byte >>= (2 * current2bits);

          int val1 = (byte & 3); // extraction des bits correspondant

          uint8_t val = val1; // pour avoir la bonne taille pour le ||

          // if (new_matrix.SNPs.size() < 10) {
          // std::cout << "val = ";
          // for (int i = 7; i >= 0; --i) {
          //     std::cout << ((val >> i) & 1);
          // }
          // std::cout << "\n";
          // }

          size_t new_byte = i / 4;
          size_t new_2bits = (i % 4) * 2;

          val <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte

          //   if (new_matrix.SNPs.size() < 10 && new_2bits) {
          //   std::cout << "val after shift = ";
          //   for (int i = 7; i >= 0; --i) {
          //       std::cout << ((val >> i) & 1);
          //   }
          //   std::cout << "\n";
          // }

          // pê un masque possible ???
          new_snp->data()[new_byte] |= val; // le || pour ne pas toucher au bits déjà set
          // if (new_matrix.SNPs.size() < 10 && i == new_nb_inds - 1) {
          // std::cout << "new byte final: ";
          // for (int i = 7; i >= 0; --i) {
          //     std::cout << ((new_snp->data()[new_byte] >> i) & 1);
          // }
          // std::cout << "\n";
          // sleep(1);
          // }
        }        
        new_matrix.SNPs.push_back(new_snp);
      }
      else {
        std::cout << "You are trying to extract SNPs vectors on disk, not implemented for now\n";
        break;
      }
      }
    return new_matrix;
  }

  /**
   * @brief A member function testing whether the SNP at index
   * is a SNPVectorDisk or a SNPVectorMemory
   *
   * @param index
   * @return true
   * @return false
   */
  bool onDisk(size_t index)
  { // TODO  check why was that written ????doesn't work, will see later
    if (std::shared_ptr<SNPVectorDisk> test = std::dynamic_pointer_cast<SNPVectorDisk>(SNPs[index]))
      return true;
    return false; //test was a null ptr, not castable to a SNPVectorDisk
  }

  const std::vector<std::shared_ptr<SNPvector>> &getSNPs() const
  {
    return SNPs;
  }

  const std::shared_ptr<SNPvector> &getSNP(size_t i) const
  {
    return SNPs[i];
  }

private:
  std::vector<std::shared_ptr<SNPvector>> SNPs;
};

#endif