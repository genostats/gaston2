#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector> // for omp reduction
#include <memory> // for shared_ptr
#include <stdexcept>  // for out of range exceptions

#include "Datastruct.h"

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
    if (nbSNPs() > 0 && nbInds() != v->nbInds())
    {
      std::cerr << "Pb loading SNP" << std::endl;
      throw std::out_of_range("Attempting to load a SNP with a different nb of individuals");
    }
    SNPs_.push_back(v);
    indStatsComputed_ = false; // si on ajoute des SNP les stats individuelles doivent être recalculées
  }

  /**
   * @brief get number of SNPs (see also nbSNPs)
   */
  unsigned int size() const { return SNPs_.size(); }
  /**
   * @brief get number of SNPs (see also size)
   */
  unsigned int nbSNPs() const { return SNPs_.size(); }
  /**
   * @brief get number of Inds
   */
  unsigned int nbInds() const {
    if(nbSNPs() == 0) {
      return 0;
    }
    return SNPs_[0]->nbInds();
  }

  // temporary func to test d°
  void deleteSNP()
  {
    SNPs_.pop_back();
  }

#pragma omp declare reduction(vec_int_plus : std::vector<int> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0))

  void compute_indStats()  {
    if(indStatsComputed_) return; // stats déjà calculées, on ne recalcule pas
                                  
    int total = SNPs_.size();
    if (total == 0)
      throw std::out_of_range("No SNPs loaded into the SNPMatrix !");

    size_t nbInds = SNPs_[0]->nbInds();
    std::vector<int> unordered_stats(nbInds * 4, 0); // DO NOT USE nbCHars ! rounded up !

#pragma omp parallel for reduction(vec_int_plus : unordered_stats)
    for (int i = 0; i < total; i++)
    { // parcourt tous les SNPs
      auto snp = SNPs_[i];
      size_t nbBytes = snp->nbChars();

      // parcourt le SNP[i], byte by byte
      for (size_t byte = 0; byte < nbBytes; byte++)
      {
        uint8_t d = snp->data()[byte];
        size_t byteoffset = byte * 4;

        // check si byte entiere à lire
        int max_ind = 4;
        if (byte == nbBytes - 1)
        {
          int reste = nbInds % 4;
          if (reste)
            max_ind = reste; // if != 0, else stays at 4 to read full byte
        }

        for (int ind = 0; ind < max_ind; ind++)
        {
          unsigned int val_plink = snp->currentMode_[0][(d >> (2 * ind)) & 3];
          unordered_stats[(byteoffset + ind) * 4 + val_plink]++;
        }
      }
    }

    std::vector<int> vecN0s;
    std::vector<int> vecN1s;
    std::vector<int> vecN2s;
    std::vector<int> vecNAs;

    for (size_t ind = 0; ind < nbInds; ind++)
    {
      auto idxnzeros = ind * 4;
      vecN0s.push_back(unordered_stats[idxnzeros]);
      vecN1s.push_back(unordered_stats[idxnzeros + 1]);
      vecN2s.push_back(unordered_stats[idxnzeros + 2]);
      vecNAs.push_back(unordered_stats[idxnzeros + 3]);
    }

    Column N0s(vecN0s);
    Column N1s(vecN1s);
    Column N2s(vecN2s);
    Column NAs(vecNAs);

    indStats_.push_back(N0s, "N0");
    indStats_.push_back(N1s, "N1");
    indStats_.push_back(N2s, "N2");
    indStats_.push_back(NAs, "NAs");

    indStatsComputed_ = true;
    // TODO : could add a checksum ?
  }

  // TODO à déplacer dans la classe DataStruct (fonction extractLines par exemple)
  DataStruct extract_indStats(const std::vector<size_t> &keep)
  {
    DataStruct filtered_stats;

    for (Column &col : indStats_.cols)
    {
      auto type = col.type();
      if (type == datatype::INT)
      {
        const auto &values = col.get<int>();
        std::vector<int> filtered;
        for (size_t idx : keep)
        {
          filtered.push_back(values->at(idx));
        }
        filtered_stats.push_back(filtered);
      }
      else if (type == datatype::FLOAT)
      {
        const auto &values = col.get<float>();
        std::vector<float> filtered;
        for (size_t idx : keep)
        {
          filtered.push_back(values->at(idx));
        }
        filtered_stats.push_back(filtered);
      }
      else if (type == datatype::DOUBLE)
      {
        const auto &values = col.get<double>();
        std::vector<double> filtered;
        for (size_t idx : keep)
        {
          filtered.push_back(values->at(idx));
        }
        filtered_stats.push_back(filtered);
      }
      else //find out how to add a custom type ? 
      {
        throw std::runtime_error("No an existing type yet");
      }
    }

    return filtered_stats;
  }

  // just a small helper f° to clrify extract_ind, not sure if usefull
  int read_ind(uint8_t byte, size_t ind_idx) {
    byte >>= (2 * (ind_idx % 4));
    return byte & 3; // extraction des bits correspondant
  }

  /**
   * @brief A member function testing whether the SNP at index
   * is a SNPvectorDisk or a SNPVectorMemory
   *
   * @param index
   * @return true
   * @return false
   */
  bool onDisk(size_t index)
  { // TODO  check why was that written ????doesn't work, will see later
    if (std::shared_ptr<SNPvectorDisk<mio::access_mode::read>> test = std::dynamic_pointer_cast<SNPvectorDisk<mio::access_mode::read>>(SNPs_[index]))
      return true;
    return false; // test was a null ptr, not castable to a SNPvectorDisk
  }

  const std::vector<std::shared_ptr<SNPvector>> &getSNPs() const
  {
    return SNPs_;
  }

  const std::shared_ptr<SNPvector> &getSNP(size_t i) const
  {
    return SNPs_[i];
  }

  const DataStruct getIndStats() const { return indStats_; }

  // TODO : see if by default possible ?
  // they need to be ordered !!!!
  void setIndStats(Column N0s, Column N1s, Column N2s, Column NAs)
  {
    indStats_.push_back(N0s);
    indStats_.push_back(N1s);
    indStats_.push_back(N2s);
    indStats_.push_back(NAs);
  }

  void setIndStats(DataStruct new_stats) {
    indStats_ = new_stats;
  }

  // compute all SNP stats
  void computeSNPStats() {
    for(auto & snp : SNPs_) {
      snp->compute_stats();
    }
  }

  // comute SNP stats for snp i with i1 <= i <= i2
  // SHOULD THIS CHANGE TO i1 <= i < i2 ?
  void computeSNPStats(size_t i1, size_t i2) {
    for(size_t i = i1; i <= i2; i++) {
      SNPs_[i]->compute_stats();
    }
  }

  
  void setMode(Mode mode) {
    for(auto & snp : SNPs_) {
      snp->setMode(mode);
    }
  }
private:
  DataStruct indStats_;
  std::vector<std::shared_ptr<SNPvector>> SNPs_;
  bool indStatsComputed_ = false;
};

#endif
