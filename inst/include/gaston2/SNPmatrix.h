#include <omp.h>

#include <fstream>    // for ifstream
#include <memory>     // for shared_ptr
#include <stdexcept>  // for out of range exceptions
#include <vector>     // for omp reduction

#include "Datastruct.h"
// TODO : temporary !!
//#include "SNPdosage.h"
//#include "SNPdosageMemory.h"
#include "SNPvector.h"

#ifndef _snpmatrix_
#define _snpmatrix_

/**
 * @brief A class keeping in a vector pointers to SNPvectors
 * (either SNPvectorDisk or SNPvectorMemory)
 *
 */
template <typename SNPvectorClass = SNPvector>
class SNPmatrix {
 public:
  /**
   * @brief function to add a shared_ptr into the vector keeping the
   * SNPs of the matrix. /!\ will check if the SNP added is matching in size
   *
   */
  void push_back(std::shared_ptr<SNPvectorClass> v) {
    // if at least one loaded, everySNPs must have same size
    if (nbSNPs() > 0 && nbInds() != v->nbInds()) {
      std::cerr << "Pb loading SNP" << std::endl;
      throw std::out_of_range("Attempting to load a SNP with a different nb of individuals");
    }
    v->setMode(mode_);
    SNPs_.push_back(v);
    indStatsComputed_ = false;  // si on ajoute des SNP les stats individuelles doivent être recalculées
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
    if (nbSNPs() == 0) {
      return 0;
    }
    return SNPs_[0]->nbInds();
  }

  // temporary func to test d°
  void deleteSNP() {
    SNPs_.pop_back();
  }

  SNPmatrix() {}  // default constructor to have a specialized one just after

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
  SNPmatrix(const SNPmatrix<SNPvectorClass> &other, intVec keep) {
    const std::vector<std::shared_ptr<SNPvectorClass>> otherSNPs = other.getSNPs();
    for (auto keep_idx : keep) {
      this->push_back(otherSNPs.at(keep_idx));  // at is supposed to do bound checking
    }
    
    // now inheriting the SNPStats_ of the SNPs specified in keep
    DataStruct original_snpStats = other.getSNPStats();
    snpStats_ = DataStruct(original_snpStats, keep);
    // just to get the fam stats in
    indStats_ = other.getIndStats();
    // Still not automatically computing indStats back,
    // so c° is setting indStatsComputed_ to false by default
  }

  // #pragma omp declare reduction(vec_int_plus : std::vector<int> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0))

  // Stocking in Datastruct indStats_ the number of occurrences of 0, 1, 2, and NAs for all individuals,
  // accross all SNPs currently in the SNPmatrix (need to be recalled if another SNP is pushed back after)
  void compute_indStats(bool force = false) {
    if (!force && indStatsComputed_) return;  // stats déjà calculées, on ne recalcule pas

    size_t nbSNPs = SNPs_.size();

    if (nbSNPs == 0)
      throw std::out_of_range("No SNPs loaded into the SNPMatrix !");

    size_t nbInds = SNPs_[0]->nbInds();
    // a vector keeping for every inds NOs, N1s, N2s and NAs
    // filled with 0 by default
    std::vector<int> unordered_stats(nbInds * 4, 0);  // DO NOT USE nbCHars ! rounded up !

    // #pragma omp parallel for reduction(vec_int_plus : unordered_stats)
    for (size_t i = 0; i < nbSNPs; i++) {  // parcourt tous les SNPs
      SNPs_[i]->compute_indStats(unordered_stats);
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

  const std::vector<std::shared_ptr<SNPvectorClass>> &getSNPs() const {
    return SNPs_;
  }

  const std::shared_ptr<SNPvectorClass> &getSNP(size_t i) const {
    return SNPs_[i];
  }

  // get the DataStruct containing individual stats
  const DataStruct &getIndStats() const {
    // TO DEBUG : will print the type of what the SNPmat has
    //if (SNPs_.size()) {
      //auto ptr = getSNP(0);
      //std::cout << "Is current SNPmatrix containing SNPdosage ? (0 = false) "
                //<< (dynamic_cast<SNPdosage*>(ptr.get()) != nullptr)
                //<< "\n";  
      //std::cout << "Is current SNPmatrix containing SNPdosageMemory ? (0 = false) "
      //<< (dynamic_cast<SNPdosageMemory*>(ptr.get()) != nullptr)
      //<< "\n";   }
      return indStats_;
  }

  // get the DataStruct containing snp stats
  const DataStruct &getSNPStats() const { return snpStats_; }

  // TODO : see if by default possible ?
  // they need to be ordered !!!!
  void setIndStats(Column N0s, Column N1s, Column N2s, Column NAs) {
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

  // for the toSNPmatrix* functions and the extractInds* f° to populate through
  // snpStats_ with at least bim content
  void setSnpStats(DataStruct new_snps) {
    snpStats_ = new_snps;
  }

  // compute all SNP stats
  void computeSNPStats() {
    for (auto &snp : SNPs_) {
      snp->compute_stats();
    }
  }

  // comute SNP stats for snp i with i1 <= i <= i2
  // SHOULD THIS CHANGE TO i1 <= i < i2 ?
  void computeSNPStats(size_t i1, size_t i2) {
    if (i1 > i2) std::swap(i1, i2);
    if (i2 >= SNPs_.size())
      throw std::out_of_range("Out of range [computeSNPStats]");
    for (size_t i = i1; i <= i2; i++) {
      SNPs_[i]->compute_stats();
    }
  }

  void setMode(Mode mode) {
    for (auto &snp : SNPs_) {
      snp->setMode(mode);
    }
    mode_ = mode;
  }

  // TODO (to think) there might be a problem if SNPs are not all in the same mode...
  // possible solution : enforce mode when push_back is done ?
  inline Mode mode() {
    return mode_;
  }

  void readFamFile(std::string famFile) {
    std::ifstream in(famFile);
    if (!in.good())
      throw std::runtime_error("Can't open fam file");
    std::vector<datatype> colTypes = {datatype::STRING, datatype::STRING, datatype::STRING, datatype::STRING, datatype::INT, datatype::INT};
    std::vector<std::string> colNames = {"famid", "id", "father", "mother", "sex", "pheno"};
    indStats_ = DataStruct(colTypes, colNames);
    indStats_.readFile(in);
  }

  void readBimFile(std::string bimFile) {
    std::ifstream in(bimFile);
    if (!in.good())
      throw std::runtime_error("Can't open bim file");
    std::vector<datatype> colTypes = {datatype::STRING, datatype::STRING, datatype::INT, datatype::DOUBLE, datatype::STRING, datatype::STRING};
    std::vector<std::string> colNames = {"chr", "id", "pos", "dist", "A1", "A2"};
    snpStats_ = DataStruct(colTypes, colNames);
    snpStats_.readFile(in);
  }

 private:
  // stats and informations
  DataStruct indStats_;  // will contain fam file + statistiques
  DataStruct snpStats_;  // will contain bim file + ?
  std::vector<std::shared_ptr<SNPvectorClass>> SNPs_;
  bool indStatsComputed_ = false;
  Mode mode_ = RAW_VALUES;
};

#endif
