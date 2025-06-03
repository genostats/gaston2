#include <cstdint>
#include <cstddef>
#include <omp.h>

#include <math.h>    // for sqrt...
#include <stdexcept> // for out of range exceptio
#include <iostream>
#include <vector>    //for contingency !!
#include <algorithm> // for omp reduction ?
#include "debug.h"

#ifndef _SNPvector_
#define _SNPvector_

static uint8_t N0[256] = {
    4, 3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2,
    3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1,
    3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1,
    3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1,
    3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

static uint8_t N1[256] = {
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1,
    2, 3, 2, 2, 3, 4, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2,
    1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1,
    1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0};

enum Mode
{
  PLINK,                 // .bed file : g = {0, 1, 3 and 2 = NA}
  CENTERED,              // g -= mu
  STANDARDIZED_MU_SIGMA, // g = (g-mu)/sd
  STANDARDIZED_P,        // g = (g - 2⁼p) / sqrt(2*p*(1-p)) with p = mu/2
  CUSTOM 
};

/**
 * @brief An abstract class instanciated through SNPvectorMemory or SNPvectorDisk
 *
 * Stores a bit vector, could be in memory or in a memory_mapped file
 *
 */
class SNPvector {

protected: // can be accessed also by class inheriting

  // nbInds and mode shared by all inheriting classes
  // (cf constructors below)

  const size_t nbInds_;
  enum Mode mode_ = Mode::PLINK;

  /* Containing N0, N1, N2, NAs on the whole SNP,
  following Plink format
  populated by compute_stats()*/

  unsigned int stats_[4] = {0, 0, 0, 0};
  bool stats_set_ = false;

  double mu_ = 0;
  double sigma_ = 0;

  // just a small helper f° to extract an individual genotype from a byte
  inline uint8_t read_ind(uint8_t byte, size_t ind_idx) {
    byte >>= (2 * (ind_idx % 4));
    return byte & 3; // extraction des bits correspondant
  }

  // an array for transformed genotype (either centered, standardized, etc)
  // values can be changed by function setMode
  // default values = "raw genotypes" with a 3 for NA (could be set to NaN ? to be considered)
  // (see also in setMode below)
  double g_trans[4] = {0, 3, 1, 2};

  // constructor allowing to set nbInds and mode (to be called by derived classes)
  SNPvector(size_t nbInds, Mode mode = Mode::PLINK) : nbInds_(nbInds), mode_(mode) {}

  inline void computeMode() { 
    switch(mode_) {
      case Mode::PLINK: {
        g_trans[0] = 0; 
        g_trans[1] = 3;   // or NAN ?
        g_trans[2] = 1;
        g_trans[3] = 2;
        break;
      }
      case Mode::CENTERED: {
        g_trans[0] = 0 - mu_; 
        g_trans[1] = 0;
        g_trans[2] = 1 - mu_;
        g_trans[3] = 2 - mu_;
        break;
      }
      case Mode::STANDARDIZED_MU_SIGMA: {
        g_trans[0] = (0 - mu_)/sigma_; 
        g_trans[1] = 0;
        g_trans[2] = (1 - mu_)/sigma_;
        g_trans[3] = (2 - mu_)/sigma_;
        break;
      }
      case Mode::STANDARDIZED_P: {
        double s = sqrt( mu_*(1 - mu_/2) ); // sqrt 2p(1-p)
        g_trans[0] = (0 - mu_)/s; 
        g_trans[1] = 0;
        g_trans[2] = (1 - mu_)/s;
        g_trans[3] = (2 - mu_)/s;
        break;
      }
      default:
        throw std::runtime_error("Private function computeMode called with bad mode value");
    }
  }

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
   * returning a pointer to the front of the vector
   *
   * @return uint8_t*
   */
  // virtual uint8_t *data() = 0;
  virtual const uint8_t *data() const = 0;

  // nombre d'individus
  size_t nbInds() const {
    return nbInds_;
  }

  // BY DEFAULT, PLINK format, with 01 = missing genotype
  // double currentMode_[5][4] = {{0, 3, 1, 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 2, 3}, {0, 0, 0, 0}};

  // To use set mode and compute g_trans value accordingly
  void setMode(Mode mode) { 
    if(mode == mode_) return; // didn't change, do nothing
    if(mode == Mode::CUSTOM) 
      throw std::runtime_error("Can't set custom mode this way");
    mode_ = mode;
    computeMode();
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

  const unsigned int *getStats() const {
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
    computeMode(); // recompute mode
  }

  void setSigma(double sigma) {
    sigma_ = sigma;
    computeMode(); // recompute mode
  }
  
  void setMuSigma(double mu, double sigma) {
    mu_ = mu;
    sigma_ = sigma;
    computeMode(); // recompute mode
  }

  // dummy function summing the n first individuals from the SNP, if no n specified, summing all SNPs available
  double sum(int n = -1) {
    double S = 0;
    const uint8_t *d = data();
    size_t nc = nbChars(); // max nb of bytes
    int cptr = 0;
    if (n == -1)
      n = nbInds(); // set to max nb of bits available

    for (size_t i = 0; (cptr < n && i < nc); i++) {
      uint8_t byte = d[i];
      for (unsigned int bits = 0; (cptr < n && bits < 4); bits++) {
        unsigned int val = (byte >> (2 * bits)) & 3; // Extract the SNP (2 bits)
        S += g_trans[val];
        cptr++;
      }
    }
    return S;
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
    if(stats_set_) return;

    size_t nbc_m1 = nbChars() - 1;

    /* FIRST : filling up stats_ with N0s, N1s, N2s, and NAs with PLINK translation*/
    stats_[0] = stats_[1] = stats_[2] = stats_[3] = 0;

    // all bytes excepted last byte
    for (size_t i = 0; i < nbc_m1; i++) {
      uint8_t d = data()[i];
      stats_[0] += N0[d];
      stats_[2] += N0[255 - d]; // raw = 3
      stats_[3] += N1[d];       // raw = 1
      stats_[1] += N1[255 - d]; // raw = 2
    }

    // last byte treated separately,
    // stopping mid-byte if necessary:
    unsigned int g[4] = {0, 3, 1, 2};
    uint8_t d = data()[nbc_m1];
    unsigned int BitsInLastByte = (nbInds() & 3)?(nbInds() & 3):4;
    while(BitsInLastByte > 0) {
      unsigned int val_plink = g[d&3];
      stats_[val_plink]++;
      BitsInLastByte--;
      d >>= 2; // 1 shift par loop
    }

    // stats are set !!
    stats_set_ = true;

    // THEN : computing mu and sigma (according to arguments)
    compute_mu_sigma(set_mu, set_sigma);

    /* FINALLY : updating the "mode" enum that acts as a filter
     to get calculated via mu_and sigma_*/
    computeMode();
  }

  // for scalar product :
  // CAVEAT: stats are supposed set! 
  template<typename scalar_t = double>
  inline scalar_t LD(const SNPvector &other) const {
     
    size_t nbi = nbInds(); 
    if (nbi != other.nbInds())
       throw std::runtime_error("Mismatch in the nb of individuals between the 2 SNPs !");
    scalar_t LD = 0;
    scalar_t gg[16];
    gg[1] = gg[4] = gg[5] = gg[6] = gg[7] = gg[9] = gg[13] = 0;

    gg[0] = (-mu_) * (-(other.mu_));
    gg[2] = (-mu_) * (1. - (other.mu_));
    gg[3] = (-mu_) * (2. - (other.mu_));

    gg[8] = (1. - mu_) * (-(other.mu_));
    gg[10] = (1. - mu_) * (1. - (other.mu_));
    gg[11] = (1. - mu_) * (2. - (other.mu_));

    gg[12] = (2. - mu_) * (-(other.mu_));
    gg[14] = (2. - mu_) * (1. - (other.mu_));
    gg[15] = (2. - mu_) * (2. - (other.mu_));

    // PB ! what if v == 0 ?
    // TODO return NaN directly
    scalar_t v = sigma_ * other.sigma_;

    auto data1 = data();
    auto data2 = other.data();
    size_t nbc_m1 = nbChars() - 1;

    for (size_t i = 0; i < nbc_m1; i++) {
      uint8_t g1 = data1[i]; // je récup les ièmes char
      uint8_t g2 = data2[i];

      for (unsigned int ss = 0; ss < 4; ss++) { 
        LD += gg[ ((g1&3)*4) + (g2&3) ];
        g1 >>= 2;
        g2 >>= 2;
      }
    }

    // gaston pouvait supposer que le dernier char était borné par des NA
    // on ne peut pas le faire ici : on traite le dernier char à part
    // (étonnament ça fait perdre plusieurs microsecondes)
    // number of bits to read on last byte
    // cette cabriole gagne un peu de temps... 
    unsigned int BitsInLastByte = (nbi & 3)?(nbi & 3):4;

    uint8_t g1 = data1[nbc_m1];
    uint8_t g2 = data2[nbc_m1];
    for (unsigned int ss = 0; ss < BitsInLastByte; ss++) {
      LD += gg[ ((g1&3)*4) + (g2&3) ];
      g1 >>= 2;
      g2 >>= 2;
    }

    scalar_t r = LD / (v * (nbi - 1)); // nbInds should be the same for both
    return r;
  }

  // function filling a 9 (unsigned) int vector corresponding to a contigency table of genotypes over 2 SNPs
  // intVec must have members size() and []
  template<typename intVec>
  void contingency(const SNPvector &other, intVec & contingencyTable) {
    if(contingencyTable.size() < 9) throw std::runtime_error("In contingency, contigencyTable is too short");
    // the "raw" contigency table, 4x4 seen as a long 16 elts array
    unsigned int table[16] = {0}; // initialisée à 0
    size_t nbc_m1 = nbChars() - 1;
    auto data1 = data();
    auto data2 = other.data();
    for (size_t i = 0; i < nbc_m1; i++) {
      uint8_t g1 = data1[i]; 
      uint8_t g2 = data2[i];
      for (unsigned int ss = 0; ss < 4; ss++) {
        table[ (g1&3)*4 + (g2&3) ]++;
        g1 >>= 2;
        g2 >>= 2;
      }
    }
    // idem LD(), dernier char traité à part
    size_t nbi = nbInds();
    unsigned int BitsInLastByte = (nbi & 3)?(nbi & 3):4;
    uint8_t g1 = data1[nbc_m1];
    uint8_t g2 = data2[nbc_m1];
    for (unsigned int ss = 0; ss < BitsInLastByte; ss++) {
      table[ (g1&3)*4 + (g2&3) ]++;
      g1 >>= 2;
      g2 >>= 2;
    }
    // on recopie les résultats en virant les 7 valeurs
    // où un des deux génotypes au moins est à NA -> table 3x3
    // linéarisée
    contingencyTable[0] = table[0];
    contingencyTable[1] = table[2];
    contingencyTable[2] = table[3];
    contingencyTable[3] = table[8];
    contingencyTable[4] = table[10];
    contingencyTable[5] = table[11];
    contingencyTable[6] = table[12];
    contingencyTable[7] = table[14];
    contingencyTable[8] = table[15];
  }

  // on va incrémenter un vecteur de scalar_t [n'importe quoi qui a un opérateur [] et size())
  // size = n * (n + 1) / 2 = la moitié d'une matrice symmétrique, diagonale incluse
  // si on voit V comme une matrice c'est V += SNP . SNP'
  template<typename scalar_t, typename vectorType> 
  void tcrossprod0(vectorType & V) {
    if(V.size()*2 != nbInds_ * (nbInds_ + 1)) throw std::runtime_error("tcrossprod, V has not the right size");

    auto DATA = data();
    size_t nbc_m1 = nbChars() - 1;
   
    // value of transformed *genotypes* 0, 1, 2, skipping NA (always 0)
    scalar_t v0 = (scalar_t) g_trans[0];
    scalar_t v1 = (scalar_t) g_trans[2];
    scalar_t v2 = (scalar_t) g_trans[3];
    scalar_t VALUES[16] = 
      { v0*v0, 0, v0*v1, v0*v2,
            0, 0,     0,     0,
        v0*v1, 0, v1*v1, v1*v2,
        v0*v2, 0, v1*v2, v2*v2 } ;
    // loop j1 to nbChar - 1
    for(size_t j1 = 0; j1 < nbc_m1; j1++) {
      // l'indice où écrire dans V
      size_t k = 2*j1 * (4*j1 + 1);
      uint8_t g1 = DATA[j1]; 
      for(unsigned int ss1 = 0; ss1 < 4; ss1++) {
        scalar_t * VALUES1 = VALUES + (g1&3) * 4;  // le début de la ligne qui correspond à (g1&3)
        g1 >>= 2;
        // on boucle j2 de 0 à j1 - 1
        for(size_t j2 = 0; j2 < j1; j2++) {
          uint8_t g2 = DATA[j2];
          for(unsigned int ss2 = 0; ss2 < 4; ss2++) {
            V[k++] += VALUES1[ g2&3 ];
            g2 >>= 2;
          }
        }
        // j2 = j1 est traité à part
        uint8_t g2 = DATA[j1];
        for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
          V[k++] += VALUES1[ g2&3 ];
          g2 >>= 2;
        }
      }
    }
    // j1 = nbChar - 1 est traité à part. C'est essentiellement la même boucle
    // (on a hardcoded j1 = nbc_m1) 
    // mais on s'arrête à BitsInLastByte pour ss1. Cette boucle ne peut être déroulée à 
    // la compilation.
    unsigned int BitsInLastByte = (nbInds_ & 3)?(nbInds_ & 3):4;
    size_t k = 2*nbc_m1 * (4*nbc_m1 + 1);
    uint8_t g1 = DATA[nbc_m1];

    for(unsigned int ss1 = 0; ss1 < BitsInLastByte; ss1++) {
      scalar_t * VALUES1 = VALUES + (g1&3) * 4;  // le début de la ligne qui correspond à (g1&3)
      g1 >>= 2;
      // on boucle j2 de 0 à j1 - 1
      for(size_t j2 = 0; j2 < nbc_m1; j2++) {
        uint8_t g2 = DATA[j2];
        for(unsigned int ss2 = 0; ss2 < 4; ss2++) {
          V[k++] += VALUES1[ g2&3 ];
          g2 >>= 2;
        }
      }
      // j2 = j1 = nbc_m1 est traité à part
      uint8_t g2 = DATA[nbc_m1];
      for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
        V[k++] += VALUES1[ g2&3 ];
        g2 >>= 2;
      }
    }
  }

  template<typename scalar_t, typename vectorType> 
  void tcrossprod(vectorType & V) {
    if(V.size()*2 != nbInds_ * (nbInds_ + 1)) throw std::runtime_error("tcrossprod, V has not the right size");

    auto DATA = data();
    size_t nbc_m1 = nbChars() - 1;
   
    // value of transformed *genotypes* 0, 1, 2, skipping NA (always 0)
    const scalar_t v0 = (scalar_t) g_trans[0];
    const scalar_t v1 = (scalar_t) g_trans[2];
    const scalar_t v2 = (scalar_t) g_trans[3];

    // loop j1 to nbChar - 1
    uint8_t previous_genotype_1 = 4; // impossible value
#pragma omp parallel for firstprivate(previous_genotype_1)
    for(size_t j1 = 0; j1 < nbc_m1; j1++) {
      // l'indice où écrire dans V : début de la ligne 4j1, on a écrit 4j1 * (4j1 + 1) / 2 éléments
      size_t k = 2*j1 * (4*j1 + 1);
      uint8_t g1 = DATA[j1]; 

      for(unsigned int ss1 = 0; ss1 < 4; ss1++) {
        uint8_t current_genotype_1 = g1&3;

        // skip if genotype is NA / all H1 = 0...
        if(current_genotype_1 == (uint8_t) 1) {
          k += 4*j1 + ss1 + 1;
          g1 >>= 2;
          continue;
        }

// TODO tester de précalculer les trois tableaux H1 possibles...
        scalar_t v_g1 = g_trans[ current_genotype_1 ];
        scalar_t H1[32];
        if(current_genotype_1 != previous_genotype_1) { // update H1 only if v_g1 did change
          const scalar_t v_g1_v0 = v_g1*v0;
          const scalar_t v_g1_v1 = v_g1*v1;
          const scalar_t v_g1_v2 = v_g1*v2;
          H1[0]  = v_g1_v0;  H1[1] = v_g1_v0;  H1[2]  = 0;        H1[3]  = v_g1_v0; 
          H1[4]  = v_g1_v1;  H1[5] = v_g1_v0;  H1[6]  = v_g1_v2;  H1[7]  = v_g1_v0; 
          H1[8]  = v_g1_v0;  H1[9] = 0;        H1[10] = 0;        H1[11] = 0; 
          H1[12] = v_g1_v1; H1[13] = 0;        H1[14] = v_g1_v2;  H1[15] = 0; 
          H1[16] = v_g1_v0; H1[17] = v_g1_v1;  H1[18] = 0;        H1[19] = v_g1_v1; 
          H1[20] = v_g1_v1; H1[21] = v_g1_v1;  H1[22] = v_g1_v2;  H1[23] = v_g1_v1; 
          H1[24] = v_g1_v0; H1[25] = v_g1_v2;  H1[26] = 0;        H1[27] = v_g1_v2; 
          H1[28] = v_g1_v1; H1[29] = v_g1_v2;  H1[30] = v_g1_v2;  H1[31] = v_g1_v2; 
        }
        // le tableau H1 est construit de façon que 
        // H1[2*x] et H1[2*x+1] sont les valeurs des g_trans[geno1] * g_trans[geno2] et g_trans[geno1] * g_trans[geno2']
        // où x = 0 à 15 (4 bits) est la concaténation de geno2 (poids faible) et geno2' (poids fort)

        g1 >>= 2;
        previous_genotype_1 = current_genotype_1;
        // current_genotype_1 = g1&3;

        // on boucle j2 de 0 à j1 - 1
        for(size_t j2 = 0; j2 < j1; j2++) {
          uint8_t g2 = DATA[j2];
          // ou coupe g2 en deux morceaux
          uint8_t gg2 = (g2&15)*2; // apparemment le compilateur n'optimise pas ça si on le met dans les crochets ci-dessous...
          V[k++] += H1[ gg2 ];
          V[k++] += H1[ gg2 + 1 ];
          g2 >>= 4;
          V[k++] += H1[ g2*2 ];
          V[k++] += H1[ g2*2 + 1 ];
        }
        // j2 = j1 est traité à part
        uint8_t g2 = DATA[j1];
        for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
          V[k++] += v_g1 * g_trans[ g2&3 ];
          g2 >>= 2;
        }
      }
    }
    // j1 = nbChar - 1 est traité à part. C'est essentiellement la même boucle
    // (on a hardcoded j1 = nbc_m1) 
    // mais on s'arrête à BitsInLastByte pour ss1. Cette boucle ne peut être déroulée à 
    // la compilation.
    unsigned int BitsInLastByte = (nbInds_ & 3)?(nbInds_ & 3):4;
    size_t k = 2*nbc_m1 * (4*nbc_m1 + 1); // idem ci-dessus avec j1 = nbc_m1
    uint8_t g1 = DATA[nbc_m1];

    for(unsigned int ss1 = 0; ss1 < BitsInLastByte; ss1++) {
      uint8_t current_genotype_1 = g1&3;

      // skip if genotype is NA / all H1 = 0...
      if(current_genotype_1 == (uint8_t) 1) {
        k += 4*nbc_m1 + ss1 + 1;
        g1 >>= 2;
        continue;
      }

      scalar_t v_g1 = g_trans[ current_genotype_1 ];
      scalar_t H1[32];
      if(current_genotype_1 != previous_genotype_1) { // update H1 only if v_g1 did change
        const scalar_t v_g1_v0 = v_g1*v0;
        const scalar_t v_g1_v1 = v_g1*v1;
        const scalar_t v_g1_v2 = v_g1*v2;
        H1[0]  = v_g1_v0;  H1[1] = v_g1_v0;  H1[2]  = 0;        H1[3]  = v_g1_v0; 
        H1[4]  = v_g1_v1;  H1[5] = v_g1_v0;  H1[6]  = v_g1_v2;  H1[7]  = v_g1_v0; 
        H1[8]  = v_g1_v0;  H1[9] = 0;        H1[10] = 0;        H1[11] = 0; 
        H1[12] = v_g1_v1; H1[13] = 0;        H1[14] = v_g1_v2;  H1[15] = 0; 
        H1[16] = v_g1_v0; H1[17] = v_g1_v1;  H1[18] = 0;        H1[19] = v_g1_v1; 
        H1[20] = v_g1_v1; H1[21] = v_g1_v1;  H1[22] = v_g1_v2;  H1[23] = v_g1_v1; 
        H1[24] = v_g1_v0; H1[25] = v_g1_v2;  H1[26] = 0;        H1[27] = v_g1_v2; 
        H1[28] = v_g1_v1; H1[29] = v_g1_v2;  H1[30] = v_g1_v2;  H1[31] = v_g1_v2; 
      }
      g1 >>= 2;
      previous_genotype_1 = current_genotype_1;

      // on boucle j2 de 0 à j1 - 1 = nbc_m1 - 1
      for(size_t j2 = 0; j2 < nbc_m1; j2++) {
        // ici on a k == ((4 * nbc_m1 + ss1) * (4 * nbc_m1 + ss1 + 1)) / 2 + 4*j2 ;
        uint8_t g2 = DATA[j2];
        // ou coupe g2 en deux morceaux
        uint8_t gg2 = (g2&15)*2; // apparemment le compilateur n'optimise pas ça si on le met dans les crochets ci-dessous...
        V[k++] += H1[ gg2 ];
        V[k++] += H1[ gg2 + 1 ];
        g2 >>= 4;
        V[k++] += H1[ g2*2 ];
        V[k++] += H1[ g2*2 + 1 ];
      }
      // j2 = j1 = nbc_m1 est traité à part
      // ici k == ((4 * nbc_m1 + ss1) * (4 * nbc_m1 + ss1 + 1)) / 2 + 4*nbc_m1;
      uint8_t g2 = DATA[nbc_m1];
      for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
        V[k++] += v_g1 * g_trans[ g2&3 ];
        g2 >>= 2;
      }
    }
  }


  class Iterator {
  private:
    size_t currentChar;     // ii dans le code RV, correspond au byte sur lequel je suis
    size_t current2bits;    // ss dans le code RV, correspond au bit (0...3) dans le byte
    SNPvector &iterated; // link to mother class
    const double * values;     // values according to current mode of iterated

  public:
    // le constructeur vide démarre à 0, et récupère l'instance qui l'appelle
    Iterator(SNPvector &it) : currentChar(0), current2bits(0), iterated(it), values(it.values()) {}

    Iterator(size_t ind, SNPvector &it) : currentChar(ind / 4), current2bits(ind % 4), iterated(it), values(it.values()) {
      if (currentChar > iterated.nbChars())
        throw std::out_of_range("End Iterator is too far !");
    }

    // TODO est-ce qu'on ne peut pas avoir le byte shifté dans la classe et faire seulement un >>= 2 à l'appel de ++
    //      (et un update tous les 4 appels) ?

    // operateur * const : renvoie la valeur 2bits par 2bits
    double operator*() {
      uint8_t byte = iterated.data()[currentChar];
      byte >>= (2 * current2bits);   // pour avoir les 2bits de poids faible
      return values[byte&3];
    }

    // Operateur ++ : passe à la valeur suivante, PRE-INCRÉMENTATION
    Iterator &operator++() {
      current2bits++;
      if (current2bits > 3) {
        current2bits = 0;
        currentChar++;
      }
      return *this;
    }

    // Opérateur de comparaison entre deux itérateurs
    bool operator!=(const Iterator &other) const { return (currentChar != other.currentChar) || (current2bits != other.current2bits); }
  };

  // begin() : renvoie un itérateur
  Iterator begin() { return Iterator(*this); }

  /* end() : renvoie un itérateur qui "pointe vers la fin"
  aka après le dernier individu*/
  Iterator end() { return Iterator(nbInds(), *this); }
};

#endif
