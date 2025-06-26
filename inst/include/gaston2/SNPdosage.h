#include <cstdint>
#include <cstddef>
#include <omp.h>

#include <math.h>    // for sqrt...
#include <stdexcept> // for out of range exceptio
#include <iostream>
#include <vector>    //for contingency !!
#include <algorithm> // for omp reduction ?

#include "mode.h"
#include "chrType.h"
#include "debug.h"

#ifndef _SNPdosage_
#define _SNPdosage_

/**
 * @brief An abstract class instanciated through SNPvectorMemory or SNPvectorDisk
 *
 * Stores a bit vector, could be in memory or in a memory_mapped file
 *
 */
class SNPdosage {

protected: // can be accessed also by class inheriting

  // nbInds and mode shared by all inheriting classes
  // (cf constructors below)

  const size_t nbInds_;
  enum Mode mode_ = Mode::RAW_VALUES;
  enum chrType chr_type_ = chrType::UNKNOWN;

    /* Containing N0, N1, N2, NAs on the whole SNP,
  following Plink format
  populated by compute_stats()*/

  int stats_[4] = {0, 0, 0, 0};
  bool stats_set_ = false;

  double mu_ = 0;
  double sigma_ = 0;

  // an array for transformed genotype (either centered, standardized, etc)
  // values can be changed by function setMode
  // default values = "raw genotypes" with a 3 for NA (could be set to NaN ? to be considered)
  // (see also in setMode below)
  double g_trans[4] = {0, 3, 1, 2};

  // constructor allowing to set nbInds and mode (to be called by derived classes)
  SNPdosage(size_t nbInds, Mode mode = Mode::RAW_VALUES) : nbInds_(nbInds), mode_(mode) {}


  // n'utilise pas le mode courant, et multiplie l'écart type par scale (typiquement sqrt(nbSNPs))
  // deux modes possibles seulement
  inline void computeScaledMode(Mode mo, double scale) { 
    switch(mo) {
      case Mode::STANDARDIZED_MU_SIGMA: {
        double s = sigma_ * scale;
        g_trans[0] = (0 - mu_)/s; 
        g_trans[1] = 0;
        g_trans[2] = (1 - mu_)/s;
        g_trans[3] = (2 - mu_)/s;
        break;
      }
      case Mode::STANDARDIZED_P: {
        double s = sqrt( mu_*(1 - mu_/2) ) * scale; // sqrt 2p(1-p)
        g_trans[0] = (0 - mu_)/s; 
        g_trans[1] = 0;
        g_trans[2] = (1 - mu_)/s;
        g_trans[3] = (2 - mu_)/s;
        break;
      }
      default:
        throw std::runtime_error("Private function computeScaledMode called with bad mode value");
    }
  }

public:
  /**
   * @brief pure virtual function,
   * returning a pointer to the front of the data (vector or ptr to file data)
   *
   * @return float*
   */
  // virtual float *data() = 0;
  virtual const float *data() const = 0;

  // nombre d'individus
  size_t nbInds() const {
    return nbInds_;
  }

  // To use set mode and compute g_trans value accordingly
  // will NOT recompute mu & sigma, need to call 
  // compute_stats or compute_mu_sigma directly if you want that
  void setMode(Mode mode) {
    if(mode == mode_) return; // didn't change, do nothing
    if(mode == Mode::CUSTOM) 
      throw std::runtime_error("Can't set custom mode this way");
    if(mode == Mode::STANDARDIZED_P)
      throw std::logic_error("You should not be trying to read a standardized_p dosage");
    mode_ = mode;
    // no need for a computed mode because we can't precalculate 4 values
    // like we did in SNPvector
  }

  // Will set the mode to CUSTOM !!
  // Standard error is *multiplied* by scale
  // This is used in GRM computations (and nowhere else I think),
  // - because it should be slightly faster to pre-scale the SNPs than to divide the
  //   matrix by sqrt(nbSNPs) or sqrt(nbSNPs - 1) at the end of the computation.
  // - because weighting SNPs is one of my long delayed projects...
  // Should we add a member 'scale' with default value = 1 ?
  void setScaledMode(Mode mode, double scale) {
    computeScaledMode(mode, scale);
    mode_ = Mode::CUSTOM;
  }

  // Personnalized mode
  // To be used when you want to give the three values 
  // (missing genotype is always set to zero)
  template<typename vector>
  void setMode(const vector & values) {
    if(values.size() != 3) throw std::runtime_error("in setMode, size should be three");
    g_trans[0] = values[0];
    g_trans[1] = 0; // missing, always zero
    g_trans[2] = values[1];
    g_trans[3] = values[2];
    mode_ = Mode::CUSTOM;
  }

  // get mode
  Mode mode() const {
    return mode_;
  }


  // -------------------- chrType -----------------------------------
  void setChrType(chrType ty) {
    chr_type_ = ty;
  }

  chrType getChrType() {
    return chr_type_;
  }

  // ----------------------------------------------------------------
  // get values
  const double * values() const {
    return g_trans;
  }

  // checked if stats are set
  bool stats_set() const {
    return stats_set_;
  }

  /**
   * @brief function calculating the size of the vector
   * from the number of individuals/samples
   *
   * @return size_t
   */
  size_t nbChars() const {
    const size_t n = nbInds();
    return n / 4 + ((n % 4 == 0u) ? 0 : 1);
  }

  const int *getStats() const {
    return stats_;
  }

  double getMu() const {
    return mu_;
  }

  double getSigma() const {
    return sigma_;
  }

  void setMu(double mu) {
    mu_ = mu;
  }

  void setSigma(double sigma) {
    sigma_ = sigma;
  }
  
  void setMuSigma(double mu, double sigma) {
    mu_ = mu;
    sigma_ = sigma;
  }

  // This function overwrites mu and/or sigma, according to the arguments
  void compute_mu_sigma(bool set_mu = true, bool set_sigma = true) {
    if(!set_mu && !set_sigma) return;    // do nothing
    if(!stats_set_) compute_stats(false, false); // compute stats first (could also throw error ?)

    double N = (double) nbInds(); // equal to ncols in gaston

    double N1s = (double) stats_[1];
    double N2s = (double) stats_[2];
    double NAs = (double) stats_[3];
    double n = N - (double) NAs;

    double m = (2 * N2s + N1s) / n;
    if(set_mu) mu_ = m;
    if(set_sigma) {
      double mu2 = m * m;
      sigma_ = std::sqrt((N1s + 4 * N2s + NAs * mu2) / (N - 1) - N / (N - 1) * mu2);
    }
  }

  // Method filling up stats[] w/ the nb of ind = 00 (...03) in the SNP.
  void compute_stats(bool set_mu = true, bool set_sigma = true) {
    // if already called, do nothing
    if(stats_set_) {
      //std::cout << "Not recomputing !\n";
      return;
    }
    // to think, most likely do a hardcall
    for (size_t ind = 0; ind < nbInds_; ind++) {
      float val = data()[ind];
      if (std::isnan(val))
        stats_[3]++; // NA == 3
      else
        // TODO : think of the sanity check 
        // no sanity check ??
        stats_[static_cast<int>(val + 0.5f)]++;
    }

    // stats are set !!
    stats_set_ = true;

    // THEN : computing mu and sigma (according to arguments)
    compute_mu_sigma(set_mu, set_sigma);
  }


  // Adding in unordered_stats 
  // (a vector with an index for N0, N1, N2, NA for every ind in the SNP)
  // the gen value for every ind for this SNP 
  // this vector has to be seen as a 'flatten' matrix with 4 rows (N0 N1 N2 NAs)
  // and nbInds columns
  void compute_indStats(std::vector<int> &unordered_stats) {  
    //#pragma omp parallel for //num_threads(4)      
    for (size_t ind = 0; ind < nbInds_; ind++) {
      float val = data()[ind];
      if (std::isnan(val))
        unordered_stats[ind * 4 + 3]++; // NA is the third value
      else {
        // https://stackoverflow.com/questions/485525/round-for-float-in-c
        int round = static_cast<int>(val + 0.5f);//std::round(val);
        if (round > -1 && round < 3) { // don't know if really usefull, and slows down
            unordered_stats[ind * 4 + round]++;
        } else { // could be NA but I prefer to throw an exception for now
            throw std::logic_error("The rounding of a dosage value failed !");
        }
      }
    }
  }

  template<typename scalar_t = double>
  inline scalar_t LD(const SNPdosage &other) const {
     throw std::logic_error("You should not call LD on a Dosage object yet !");
  }

  // function filling a 9 (unsigned) int vector corresponding to a contigency table of genotypes over 2 SNPs
  // intVec must have members size() and []
  template<typename intVec>
  void contingency(const SNPdosage &other, intVec & contingencyTable) {
    throw std::logic_error("You should not call contingency on a Dosage object yet !");
    }

  // on va incrémenter un vecteur de scalar_t [n'importe quoi qui a un opérateur [] et size())
  // size = n * (n + 1) / 2 = la moitié d'une matrice symmétrique, diagonale incluse
  // si on voit V comme une matrice c'est V += SNP . SNP'
  template<typename scalar_t, typename vectorType> 
  void tcrossprod(vectorType & V) {
    throw std::logic_error("You should not call tcrossprod on a Dosage object yet !");
  }


  class Iterator {
  private:
    size_t currentInd;     // index de l'élement actuel
    SNPdosage &iterated; // link to mother class

  public:
    // le constructeur vide démarre à 0, et récupère l'instance qui l'appelle
    Iterator(SNPdosage &it) : currentInd(0), iterated(it) {}

    Iterator(size_t ind, SNPdosage &it) : currentInd(ind), iterated(it) {
      if (currentInd > iterated.nbInds())
        throw std::out_of_range("End Iterator is too far !");
    }

    // operateur * const : renvoie la valeur individus par individus
    double operator*() {
      float val= iterated.data()[currentInd];
      switch(iterated.mode_) {
        case Mode::RAW_VALUES: {
          return val;
        }
        case Mode::CENTERED: {
          if (isnan(val)) return 0; // TODO : check, apparently there is a bit pattern, maybe it is faster ?
          return val - iterated.mu_;
        }
        case Mode::STANDARDIZED_MU_SIGMA: {
          if (isnan(val)) return 0;
          return (val - iterated.mu_) / iterated.sigma_;
        }
        case Mode::STANDARDIZED_P: {
          throw std::logic_error("You should not be trying to read a standardized_p dosage");
        }
        case Mode::CUSTOM : {
          throw std::logic_error("You should not be trying to read a custom dosage ... (for now :)");
        }
        default :
          throw std::logic_error("I don't know what you tried to do but it failed");
      }
    }

    // Operateur ++ : passe à la valeur suivante, PRE-INCRÉMENTATION
    Iterator &operator++() {
      currentInd++;
      return *this;
    }

    // Opérateur de comparaison entre deux itérateurs
    bool operator!=(const Iterator &other) const { return (currentInd != other.currentInd); }
  };

  // begin() : renvoie un itérateur
  Iterator begin() { return Iterator(*this); }

  /* end() : renvoie un itérateur qui "pointe vers la fin"
  aka après le dernier individu*/
  Iterator end() { return Iterator(nbInds(), *this); }
};

#endif
