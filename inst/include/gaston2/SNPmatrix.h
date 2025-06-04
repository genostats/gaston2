#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector>     // for omp reduction
#include <memory>     // for shared_ptr
#include <stdexcept>  // for out of range exceptions
#include <fstream>    // for ifstream
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
  size_t size() const { 
    return SNPs_.size(); 
  }
  /**
   * @brief get number of SNPs (see also size)
   */
  size_t nbSNPs() const { 
    return SNPs_.size(); 
  }
  /**
   * @brief get number of Inds
   */
  size_t nbInds() const {
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

  SNPmatrix() {} // default constructor to have a specialized one just after

  /**
   * @brief Constructor "by copy", 
   * but only copying the SNPs from SNPmatrix "other"
   * specified in "keep".
   * @param other
   * SNPmatrix acting as a reference to extract from
   * @param keep, 
   * a templated vector of index, corresponding to the SNPs 
   * we want to keep from the "other" matrix
   * @return the new SNPmatrix
   */
  template <typename intVec>
  SNPmatrix(const SNPmatrix &other, intVec keep) {
    const std::vector<std::shared_ptr<SNPvector>> otherSNPs = other.getSNPs();
    for (auto keep_idx : keep) {
      this->push_back(otherSNPs.at(keep_idx)); // at is supposed to do bound checking
    }
    //Not automatically computing indStats back,
    // so c° is setting indStatsComputed_ to false by default
  }

#pragma omp declare reduction(vec_int_plus : std::vector<int> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0))

  void compute_indStats(bool force = false)  {
    if(!force && indStatsComputed_) return; // stats déjà calculées, on ne recalcule pas

    size_t nbSNPs = SNPs_.size();
 
    // TODO is this really needed?
    if (nbSNPs == 0)
      throw std::out_of_range("No SNPs loaded into the SNPMatrix !");

    size_t nbInds = SNPs_[0]->nbInds();
    std::vector<int> unordered_stats(nbInds * 4, 0); // DO NOT USE nbCHars ! rounded up !

    const unsigned int g[4] = {0, 3, 1, 2};
#pragma omp parallel for reduction(vec_int_plus : unordered_stats)
    for (size_t i = 0; i < nbSNPs; i++) { // parcourt tous les SNPs
      auto snp = SNPs_[i];
      size_t nbBytes = snp->nbChars();

      // parcourt le SNP[i], byte by byte
      size_t nbc_m1 = snp->nbChars() - 1;
      for (size_t byte = 0; byte < nbc_m1; byte++) {
        uint8_t d = snp->data()[byte];
        size_t byteoffset = byte * 4;

        for (int ind = 0; ind < 4; ind++) {
          unsigned int val_plink = g[ d&3 ];
          unordered_stats[(byteoffset + ind) * 4 + val_plink]++;
          d >>= 2;
        }
      }
      // last byte read separately
      unsigned int BitsInLastByte = (nbInds & 3)?(nbInds & 3):4;
      uint8_t d = snp->data()[nbc_m1];
      size_t byteoffset = nbc_m1 * 4;
      for (int ind = 0; ind <  BitsInLastByte; ind++) {
        unsigned int val_plink = g[ d&3 ];
        unordered_stats[(byteoffset + ind) * 4 + val_plink]++;
        d >>= 2;
      }
    }

    // isolate columns from unordered_stats
    std::vector<int> vecN0s;
    std::vector<int> vecN1s;
    std::vector<int> vecN2s;
    std::vector<int> vecNAs;

    for (size_t ind = 0; ind < nbInds; ind++) {
      auto idxnzeros = ind * 4;
      vecN0s.push_back(unordered_stats[idxnzeros]);
      vecN1s.push_back(unordered_stats[idxnzeros + 1]);
      vecN2s.push_back(unordered_stats[idxnzeros + 2]);
      vecNAs.push_back(unordered_stats[idxnzeros + 3]);
    }

    indStats_.setColumn(Column(vecN0s), "N0");
    indStats_.setColumn(Column(vecN1s), "N1");
    indStats_.setColumn(Column(vecN2s), "N2");
    indStats_.setColumn(Column(vecNAs), "NAs");

    indStatsComputed_ = true;
  }

  const std::vector<std::shared_ptr<SNPvector>> &getSNPs() const
  {
    return SNPs_;
  }

  const std::shared_ptr<SNPvector> &getSNP(size_t i) const
  {
    return SNPs_[i];
  }

  // get the DataStruct containing individual stats
  const DataStruct & getIndStats() const { return indStats_; }

  // get the DataStruct containing snp stats
  const DataStruct & getSNPStats() const { return snpStats_; }


  // TODO : see if by default possible ?
  // they need to be ordered !!!!
  void setIndStats(Column N0s, Column N1s, Column N2s, Column NAs)
  {
    indStats_.push_back(N0s);
    indStats_.push_back(N1s);
    indStats_.push_back(N2s);
    indStats_.push_back(NAs);
    indStatsComputed_ = true;
  }

  void setIndStats(DataStruct new_stats) {
    indStats_ = new_stats;
    indStatsComputed_ = true;
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
    if(i1 > i2) std::swap(i1, i2);
    if(i2 >= SNPs_.size()) 
      throw std::out_of_range("Out of range [computeSNPStats]");
    for(size_t i = i1; i <= i2; i++) {
      SNPs_[i]->compute_stats();
    }
  }
  
  void setMode(Mode mode) {
    for(auto & snp : SNPs_) {
      snp->setMode(mode);
    }
    mode_ = mode;
  }

  // TODO (to thinkà there might be a problem if SNPs are not all in the same mode...
  // possible solution : enforce mode when push_back is done ?
  inline Mode mode() {
    return mode_;
  }

  void readFamFile(std::string famFile) {
    std::ifstream in(famFile);
    if(!in.good()) 
      throw std::runtime_error("Can't open fam file");
    std::vector<datatype> colTypes = { datatype::STRING, datatype::STRING, datatype::STRING, datatype::STRING, datatype::INT, datatype::INT };
    std::vector<std::string> colNames = { "famid", "id", "father", "mother", "sex", "pheno" };
    indStats_ = DataStruct(colTypes, colNames);
    indStats_.readFile(in);
  }

  void readBimFile(std::string bimFile) {
    std::ifstream in(bimFile);
    if(!in.good()) 
      throw std::runtime_error("Can't open bim file");
    std::vector<datatype> colTypes = { datatype::STRING, datatype::STRING, datatype::INT, datatype::DOUBLE, datatype::STRING, datatype::STRING };
    std::vector<std::string> colNames = { "chr", "id", "pos", "dist", "A1", "A2" };
    snpStats_ = DataStruct(colTypes, colNames);
    snpStats_.readFile(in);
  }

private:
  // stats and informations
  DataStruct indStats_;   // will contain fam file + statistiques
  DataStruct snpStats_;   // will contain bim file + ?
  std::vector<std::shared_ptr<SNPvector>> SNPs_;
  bool indStatsComputed_ = false;
  Mode mode_ = RAW_VALUES;
};

#endif
