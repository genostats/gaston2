#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector>
#include <memory>
#include <array>     // for compute stats function
#include <stdexcept> // for out of range exceptions

#ifndef _snpmatrix_
#define _snpmatrix_

/**
 * @brief A class keeping in a vector pointers to SNPvectors
 * (either SNPvectorDisk or SNPvectorMemory)
 *
 * TODO :and not a very good idea to leave SNPs public
 * --> will change
 */
class SNPmatrix
{
public:
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

  // générer un tableaux avec les N0s N1s N2s pour tous les individus
  // need to keep this array as member ?
  // void compute_stats()

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
    for (int i = 0; i < total; i++)
    { // parcourt tous les SNPs
      auto vec = SNPs[i];
      int j = 0;
      for (auto pa = vec->begin(); pa != vec->end(); ++pa)
      { // parcourt le SNP[i]
        if (j == Indnb)
        {
          stats[(int)*pa]++;
        }
        j++;
      }
    }
    return stats;
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
    return false;
  }

  std::vector<std::shared_ptr<SNPvector>> SNPs;
};

#endif