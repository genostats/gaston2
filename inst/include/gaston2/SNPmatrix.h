#include <omp.h>
#include <csignal>

#include <fstream>    // for ifstream
#include <memory>     // for shared_ptr
#include <stdexcept>  // for out of range exceptions
#include <vector>     // for omp reduction

#include "DataStruct.h"
#include "SNPvector.h"

#include "debug_flags.h"
#include "debug.h"

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
   */
  void push_back(std::shared_ptr<SNPvectorClass> v) {
    // if at least one loaded, everySNPs must have same size
    if (nbSNPs() > 0 && nbInds() != v->nbInds()) {
      std::cerr << "Pb loading SNP" << std::endl;
      throw std::out_of_range("Attempting to load a SNP with a different nb of individuals");
    }
    v->setMode(mode_);
    SNPs_.push_back(std::move(v));
    indStatsComputed_ = false;  // si on ajoute des SNP les stats individuelles doivent être recalculées
    snpStatsExported_ = false;  // idem, il va leur manquer les stats de ce nv SNP
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

  const std::vector<std::shared_ptr<SNPvectorClass>> & getSNPs() const {
    return SNPs_;
  }

  const std::shared_ptr<SNPvectorClass> & getSNP(size_t i) const {
    return SNPs_[i];
  }


  SNPmatrix() {}  // default constructor to have a specialized one just after

  /**
   * @brief Constructor "by copy",
   * but only copying the SNPs from SNPmatrix "other"
   * specified in "keep".
   * /!\ If the mode is changed in the new matrix,
   * it will be changed also for the snps of the old matrix
   * @param other
   * SNPmatrix acting as a reference to extract from
   * @param keep,
   * a templated vector of index, corresponding to the SNPs
   * we want to keep from the "other" matrix
   * @return the new SNPmatrix
   */
  template <typename intVec>
  SNPmatrix(const SNPmatrix<SNPvectorClass> &other, intVec keep) {
    const std::vector<std::shared_ptr<SNPvectorClass>> & otherSNPs = other.getSNPs();
    for (auto keep_idx : keep) {
      this->push_back(otherSNPs.at(keep_idx));  // at is supposed to do bound checking
    }

    // now inheriting the SNPStats_ of the SNPs specified in keep
    DataStruct original_snpStats = other.getSNPStats();
    snpStats_ = DataStruct(original_snpStats, keep);
    snpStatsExported_ = other.snpStatscomplete(); // true ONLY if other's have N0etc

    // just to get the fam stats in
    indStats_ = other.getIndStats();
    // Still not automatically computing indStats back,
    // so c° is setting indStatsComputed_ to false by default

    // propagating the mode, BEWARE will change
    mode_ = other.getMode();
  }


  /**
   * @brief Constructor concatenating 2 SNPmatrix,
   * and also their stats
   * (doesn't matter if one is on Disk and the other in Memory)
   * BUT it needs to be (TODO finish this sentence please Juliette ???)
   * @param first
   * SNPmatrix to append to
   * @param second
   * SNPmatrix that will be appended to "first"
   * @return the new concatenated SNPmatrix
   */
  SNPmatrix(const SNPmatrix<SNPvectorClass> &first, const SNPmatrix<SNPvectorClass> &second) {
    // I prefer to throw an error here, but push_back would also do it (to think...)
    if (first.nbInds() != second.nbInds())
      throw std::logic_error("You should not be concatenating 2 SNPmatrix with a different number of individuals !");
    const std::vector<std::shared_ptr<SNPvectorClass>> & firstSNPs = first.getSNPs();
    for (auto first_snp : firstSNPs) {
      this->push_back(first_snp);
    }
    const std::vector<std::shared_ptr<SNPvectorClass>> & scdSNPs = second.getSNPs();

    for (auto scd_snp : scdSNPs) {
      this->push_back(scd_snp);
    }

    // TODO : there is something smarter to do here
    // we could add the N0etc from first and second (if they both exist)
    // And have indStatsComputed_ to true
    // BUT maybe stats are never used, so not usefull

    // Because all individual stats are changed (N0, N1, N2...)
    indStatsComputed_ = false;
    // but the fam stay the same !
    indStats_ = first.indStats_; // there should be the same inds in first and scd !!!

    // For snpStats_, don't have to touch, so only appending
    snpStats_ = DataStruct(first.snpStats_, second.snpStats_);
    // N0etc Columns need to exist in both for them to be added by DataStruct c° 
    if (first.snpStatscomplete() && second.snpStatscomplete()) {
      snpStatsExported_ = true;
      #if DEBUG_CBIND
      std::cout << "When concat 2 matrices, both had complete snpStats\n";
      #endif
    }
    else {
      snpStatsExported_ = false; // the N0etc Columns are not computed nor exported
      #if DEBUG_CBIND
      std::cout << "When concat 2 matrices, snpStats were incomplete\n";
      #endif
    }
    // pas besoin de propag le chr_type vu que stocké ds chaque snp
  }

  /******************** IND STATS ***********************/
  
  /*
    #pragma omp declare reduction(vec_int_plus : std::vector<int> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0))
  */

 private:
  // helper for the next function
  template<typename scalar_t = float>
  void compute_indStats_chrType(chrType ct, std::string suffix) {
    size_t nbInds = SNPs_[0]->nbInds();
    size_t nbSNPs = SNPs_.size();

    scalar_t nb_snps_of_type = (scalar_t) nbSnpType[ (size_t) ct ];
    // a vector keeping for every inds NOs, N1s, N2s and NAs
    // filled with 0 by default
    std::vector<unsigned int> unordered_stats(nbInds * 4, 0); 

    // # pragma omp parallel for reduction(vec_int_plus : unordered_stats)
    for (size_t i = 0; i < nbSNPs; i++) {  // parcourt tous les SNPs
      if(SNPs_[i]->getChrType() == ct) SNPs_[i]->compute_indStats(unordered_stats);
    }

    // isolate columns from unordered_stats
    std::vector<int> N0s;
    std::vector<int> N1s;
    std::vector<int> N2s;
    std::vector<int> NAs;
    // + compute callrate and heterozigosity
    std::vector<scalar_t> callrate;
    std::vector<scalar_t> hz;

    for (size_t ind = 0; ind < nbInds; ind++) {
      auto idxnzeros = ind * 4;
      unsigned int n0 = unordered_stats[idxnzeros];
      unsigned int n1 = unordered_stats[idxnzeros + 1];
      unsigned int n2 = unordered_stats[idxnzeros + 2];
      unsigned int na = unordered_stats[idxnzeros + 3];
      N0s.push_back(n0);
      N1s.push_back(n1);
      N2s.push_back(n2);
      NAs.push_back(na);
      callrate.push_back( 1.0 - (scalar_t) na / nb_snps_of_type );
      hz.push_back( (scalar_t) n1 / (nb_snps_of_type - (scalar_t) na) );
    }

    // using setColumn will check if the Column exists 
    // and replace it, else it will push_back
    indStats_.setColumn(Column(N0s), "N0" + suffix);
    indStats_.setColumn(Column(N1s), "N1" + suffix);
    indStats_.setColumn(Column(N2s), "N2" + suffix);
    indStats_.setColumn(Column(NAs), "NAs" + suffix);
    indStats_.setColumn(Column(callrate), "callrate" + suffix);
    indStats_.setColumn(Column(hz), "hz" + suffix);
  }
 
 public:
 
  // populating indStats 
  template<typename scalar_t = float>
  void compute_indStats(bool force = false) {
    if (!force && indStatsComputed_) { 
      return;  // stats déjà calculées, on ne recalcule pas
    }

    size_t nbSNPs = SNPs_.size();
    if (nbSNPs == 0)
      throw std::out_of_range("No SNPs loaded into the SNPMatrix !");

    if(nbSnpType[ (size_t) chrType::AUTOSOME ] > 0) compute_indStats_chrType<scalar_t>(chrType::AUTOSOME, "");
    if(nbSnpType[ (size_t) chrType::X ] > 0) compute_indStats_chrType<scalar_t>(chrType::X, ".x");
    if(nbSnpType[ (size_t) chrType::Y ] > 0) compute_indStats_chrType<scalar_t>(chrType::Y, ".y");
    if(nbSnpType[ (size_t) chrType::MT ] > 0) compute_indStats_chrType<scalar_t>(chrType::MT, ".mt");

    indStatsComputed_ = true;
  }

  // get the DataStruct containing individual stats
  const DataStruct &getIndStats() const {
    // TO DEBUG : will print the type of what the SNPmat has
    // if (SNPs_.size()) {
    // auto ptr = getSNP(0);
    // std::cout << "Is current SNPmatrix containing SNPdosage ? (0 = false) "
    //<< (dynamic_cast<SNPdosage*>(ptr.get()) != nullptr)
    //<< "\n";
    // std::cout << "Is current SNPmatrix containing SNPdosageMemory ? (0 = false) "
    //<< (dynamic_cast<SNPdosageMemory*>(ptr.get()) != nullptr)
    //<< "\n";   }
    return indStats_;
  }

  // get the DataStruct containing individual stats
  // non const version to send to manipulate in R
  DataStruct & getIndStats() { return indStats_; }

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

  bool indStatscomplete() const {
    return indStatsComputed_;
  }

  void setindStatscomplete(bool complete) {
    indStatsComputed_ = complete;
  }

  /******************* SNP STATS ***********************/

  // get theDataStruct containing snp stats, possibly without N0s...
  // checkout exportSNPStats if thats what you want
  const DataStruct & getSNPStats() const { return snpStats_; }

  // non const version to send to manipulate in R
  DataStruct & getSNPStats() { return snpStats_; }

  // for the toSNPmatrix* functions and the extractInds* f° to populate through
  // snpStats_ with at least bim content
  void setSnpStats(DataStruct new_snps) {
    snpStats_ = new_snps;
    snpStatsExported_ = true;
  }

  // compute all SNP stats
  void computeSNPStats() {
    for (auto & snp : SNPs_) {
      snp->compute_stats();
    }
  }

  // compute SNP stats for snp i with i1 <= i <= i2
  // SHOULD THIS CHANGE TO i1 <= i < i2 ?
  void computeSNPStats(size_t i1, size_t i2) {
    if (i1 > i2) std::swap(i1, i2);
    if (i2 >= SNPs_.size())
      throw std::out_of_range("Out of range [computeSNPStats]");
    for (size_t i = i1; i <= i2; i++) {
      SNPs_[i]->compute_stats();
    }
  }

  // Adds N0s... Columns in the snpStats_ DataStruct
  template<typename scalar_t = float>
  void exportSNPStats(bool force) {
    // for values not computed (women only stats for X/Y chr)
    constexpr int na_value = -2147483648; // I can't find of a better solution for the moment

    size_t nb_inds = nbInds();

    // Add a check on if snpStats were computed already, 
    // then returns if no force
    if (!(force) && (snpStatsExported_)){
      return;
    }

    // for all SNPs get N0 N1 N2 on all individuals
    std::vector<int> N0s;
    std::vector<int> N1s;
    std::vector<int> N2s;
    std::vector<int> NAs;

    for (auto & snp : SNPs_) {
      if (snp->stats_set() == 0) {
        snp->compute_stats(); // returns if stats already up to date
      }
      const int * stats = snp->getStats();

      N0s.push_back(stats[0]);
      N1s.push_back(stats[1]);
      N2s.push_back(stats[2]);
      NAs.push_back(stats[3]);
    }
    snpStats_.setColumn(Column(N0s), "N0");
    snpStats_.setColumn(Column(N1s), "N1");
    snpStats_.setColumn(Column(N2s), "N2");
    snpStats_.setColumn(Column(NAs), "NAs");

    // extra computation for X/Y chr (if any)
    unsigned int nb_male = 0;

    std::vector<int> N0f;
    std::vector<int> N1f;
    std::vector<int> N2f;
    std::vector<int> NAf;

    if(nbSnpType[ (size_t) chrType::X ] > 0 || nbSnpType[ (size_t) chrType::Y ] > 0) {
      // vector of sex
      const std::vector<int> & sex = *indStats_.getColumn("sex").get<int>();
      if(sex.size() != nb_inds) throw std::runtime_error("Sex vector size doesn't match the number of individuals");
      // create mask to remove men
      size_t nbc = SNPs_[0]->nbChars();
      std::vector<uint8_t> mask(nbc);
      size_t jj = 0;
      for(size_t j = 0; j < nbc - 1; j++) {
        // men are 1
        if(sex[jj++] == 1) { nb_male++; mask[j] |= 3;   } // jj = 4*j
        if(sex[jj++] == 1) { nb_male++; mask[j] |= 12;  } // jj = 4*j + 1
        if(sex[jj++] == 1) { nb_male++; mask[j] |= 48;  } // jj = 4*j + 2
        if(sex[jj++] == 1) { nb_male++; mask[j] |= 192; } // jj = 4*j + 3
      }
      // last byte
      if(jj < sex.size() && sex[jj++] == 1) { nb_male++; mask[nbc - 1] |= 3;   }
      if(jj < sex.size() && sex[jj++] == 1) { nb_male++; mask[nbc - 1] |= 12;  }
      if(jj < sex.size() && sex[jj++] == 1) { nb_male++; mask[nbc - 1] |= 48;  }
      if(jj < sex.size() && sex[jj++] == 1) { nb_male++; mask[nbc - 1] |= 192; }
      
      // run throuh SNPs and compute women only stats only for X/Y chr
      for (auto & snp : SNPs_) {
        int count[4] = {0, 0, 0, 0};
        if(snp->getChrType() == chrType::X || snp->getChrType() == chrType::Y) {
          snp->compute_stats_mask(count, mask, nb_male);
          N0f.push_back(count[0]);
          N1f.push_back(count[1]);
          N2f.push_back(count[2]);
          NAf.push_back(count[3]);
        } else {
          N0f.push_back(na_value);
          N1f.push_back(na_value);
          N2f.push_back(na_value);
          NAf.push_back(na_value);
        }
      }

      // set values
      snpStats_.setColumn(Column(N0f), "N0.f");
      snpStats_.setColumn(Column(N1f), "N1.f");
      snpStats_.setColumn(Column(N2f), "N2.f");
      snpStats_.setColumn(Column(NAf), "NAs.f");
    } // ------------------------------------------ fin calcul chr X/Y
    unsigned int nb_female = nb_inds - nb_male; // we assume no other sex than 1/2 in the data.

    // now compute callrate, p, maf, hz, [hwe ?]
    std::vector<scalar_t> callrate;
    std::vector<scalar_t> p;
    std::vector<scalar_t> maf;
    std::vector<scalar_t> hz;

    size_t i = 0;
    for (auto & snp : SNPs_) {
      if(snp->getChrType() == chrType::Y) {
        // compute callrate on men only
        callrate.push_back( 1.0 - (scalar_t) (NAs[i] - NAf[i]) / (scalar_t) nb_male );
      } else {
        // compute callrate on all individuals
        callrate.push_back( 1.0 - (scalar_t) NAs[i] / (scalar_t) nb_inds );
      }

      scalar_t p_, hz_;
      if(snp->getChrType() == chrType::X) {
        // compute p assuming genetic values are 0/2 or 0/1 for men on chr X
        // that is 0 vs (1 or 2).
        // we compute 1 - freq(A1) 
        // denominator is (2N0f + N1f) + N0m where N0m = N0s - N0f, which simplifies to N0f + N1f + N0s
        // numerator is 2*(nb_female - NAf) + (nb_male - NAm) where NAm = NAs - NAf
        //           which simplifies to 2*nb_female - NAf + nb_male - NAs
        p_ = 1.0 - (scalar_t) (N0f[i] + N1f[i] + N0s[i]) / (scalar_t) (2*nb_female - NAf[i] + nb_male - NAs[i]) ;
        // compute hz on women only
        hz_ = (scalar_t) N1f[i] / (scalar_t) (nb_female - NAf[i]);
      } else {
        // computation of p and hz is straightforward
        // TODO take p from the SNPvector ?
        p_ = (scalar_t) (N1s[i] + 2*N2s[i]) / (scalar_t) (2*(nb_inds - NAs[i]));
        hz_ = (scalar_t) N1s[i] / (scalar_t) (nb_inds - NAs[i]);
      }

      p.push_back(p_);
      maf.push_back( (p_ < 0.5)?p_:(1.0 - p_) );
      hz.push_back(hz_);

      i++;
    }

    snpStats_.setColumn(Column(p), "p");
    snpStats_.setColumn(Column(maf), "maf");
    snpStats_.setColumn(Column(callrate), "callrate");
    snpStats_.setColumn(Column(hz), "hz");

    snpStatsExported_ = true;
  }

  // To know if snpStats contain N0, N1 ...
  // or were manually set
  bool snpStatscomplete() const {
    return snpStatsExported_;
  }

  // used when setSnpStats called to add stats
  // but still want the N0etc Columns for example
  void setSnpStatsComplete(bool complete) {
    snpStatsExported_ = complete;
  }

  /****************************************************************/

  // TODO (to think) there might be a problem if SNPs are not all in the same mode...
  // possible solution : enforce mode when push_back is done ?
  inline Mode mode() {
    return mode_;
  }

  Mode getMode() const {
    return mode_;
  }

  Mode getMode() {
    return mode_;
  }


  void setMode(Mode mode) {
    for (auto &snp : SNPs_) {
      snp->setMode(mode);
    }
    mode_ = mode;
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
    std::vector<datatype> colTypes = { datatype::INT, datatype::STRING, datatype::DOUBLE, datatype::INT, datatype::STRING, datatype::STRING };
    std::vector<std::string> colNames = { "chr", "id", "dist", "pos", "A1", "A2" };
    snpStats_ = DataStruct(colTypes, colNames);
    snpStats_.readFile(in);
  }

  // to call once bim file and SNPs are loaded, to set chromosome type in SNPs
  // nothing done yet for loading haploptypes...
  void setChrType() {
    // get 'chr' column
    if(!snpStats_.hasColumn("chr")) 
      throw std::runtime_error("No column 'chr' in snp stats. Load bim file before calling setChrType.");
    std::vector<int> & chr = *snpStats_.getColumn("chr").get<int>();
    size_t n = nbSNPs();
    if(n != chr.size())
      throw std::runtime_error("chr data size doesn't match SNPmatrix size");

    // wipe nbSnpType
    for(unsigned int & i : nbSnpType) i = 0;

    // constexpr size_t nbt = (size_t) chrType::UNKNOWN + 1:

#pragma omp parallel for reduction(+:nbSnpType[:std::size(nbSnpType)])
    for(size_t i = 0; i < n; i++) {
      chrType ct = intToChrType(chr[i]);
      SNPs_[i]->setChrType(ct);
      nbSnpType[ (size_t) ct ]++;
    }
  }

private:
  // stats and informations
  DataStruct indStats_;  // will contain fam file + statistics of Inds
  DataStruct snpStats_;  // will contain bim file + statistics of SNP
  std::vector<std::shared_ptr<SNPvectorClass>> SNPs_;
  bool indStatsComputed_ = false;
  bool snpStatsExported_ = false;
  Mode mode_ = RAW_VALUES;

  // number of SNP by types
  unsigned int nbSnpType[ (size_t) chrType::UNKNOWN + 1 ] = {0};
};

#endif
