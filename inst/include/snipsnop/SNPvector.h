#include <cstdint>
#include <cstddef>
#include <omp.h>

#include <math.h>       // for sqrt...
#include <stdexcept>      // for out of range exceptio
#include <iostream>
#include <vector> //for contingency !!
#include <algorithm> // for omp reduction ? 

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

  enum Mode { 
    PLINK = 0, // .bed file : (g = {0, 1, 3 and 2 = NA})
    CENTERED = 1, // (g -= mu)
    STANDARDIZED_MU_SIGMA = 3, // (g = (g-mu)/sd)
    NUMERIC = 4 // (g = {0, 1, 2 and 3 = NA}) (almost never used)
  };

/**
 * @brief An abstract class instanciated through SNPvectorMemory or SNPvectorDisk
 * 
 * Stores a bit vector, could be in memory or in a memory_mapped file
 * 
 */
class SNPvector {

  protected : // can be accessed also by class inheriting
  /* Containing N0, N1, N2, NAs on the whole SNP, 
  following Plink format
  populated by compute_stats()*/
  unsigned int stats_[4] = {0, 0, 0, 0};

  double mu_ = 0;
  double sigma_ = 0;

  public:
  /**
   * @brief pure virtual function, 
   * returning a pointer to the front of the vector
   * 
   * @return uint8_t* 
   */
  virtual uint8_t * data() = 0;
  virtual const uint8_t * data() const = 0;
  // nombre d'individus
  virtual size_t nbInds() = 0;

  // BY DEFAULT, PLINK format, with 01 = missing genotype
  double currentMode_[5][4] = {{0, 3, 1, 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 2, 3}};

  virtual void setMode(Mode mode) = 0;
  virtual double * mode() = 0;
  virtual double mode(unsigned int n) = 0;

  /**
   * @brief function calculating the size of the vector
   * from the number of individuals/samples 
   * 
   * @return size_t 
   */
  size_t nbChars() {
    const size_t n = nbInds();
    return n/4 + ((n%4 == 0u)?0:1);
  }

  unsigned int * getStats() {
    return stats_;
  }

  double getMu() {
    return mu_;
  }

  void setMu(double mu) {
    mu_ = mu;
    compute_mode();
  }

  double getSigma() {
    return sigma_;
  }

  void setSigma(double sigma) {
    sigma_ = sigma;
    compute_mode();
  }

  // dummy function summing the n first individuals from the SNP, if no n specified, summing all SNPs available
  double sum(int n = -1) {
    double S = 0;
    uint8_t * d = data();
    size_t nc = nbChars(); //max nb of bytes
    int cptr = 0;
    if (n == -1) n = nbInds(); //set to max nb of bits available
    
    for(size_t i = 0; (cptr < n && i < nc); i++) {
      uint8_t byte = d[i];
      for (unsigned int bits = 0; (cptr < n && bits < 4); bits++) {
          unsigned int val = (byte >> (2 * bits)) & 3; // Extract the SNP (2 bits)
          S += mode(val);
          cptr++;
      }
    }
    return S;
  }

// CAREFULL ! This function WILL overwrite mu and sigma !
  void compute_mu_sigma() {
    double N = nbInds(); //equal to ncols in gaston

    unsigned int N1s = stats_[1];
    unsigned int N2s = stats_[2];
    unsigned int NAs = stats_[3];
    double n = N - NAs;

    mu_ = (2 * N2s + N1s) / n;
    double mu2 = mu_ * mu_;
    sigma_ = sqrt(( N1s + 4 * N2s + NAs * mu2) / (N - 1) - N / (N - 1) * mu2);
  }

  void compute_mode() {

     // Centered mode Plinked
     currentMode_[1][0] = -mu_;
     currentMode_[1][1] = 0;
     currentMode_[1][2] = 1 - mu_;
     currentMode_[1][3] = 2 - mu_;

     // Centré réduit
     currentMode_[3][0] = currentMode_[1][0] / sigma_;
     currentMode_[3][1] = 0;
     currentMode_[3][2] = currentMode_[1][2] / sigma_;
     currentMode_[3][3] = currentMode_[1][3] / sigma_;
  }


  // TODO : redo this funciton
  // Method filling up stats[] w/ the nb of ind = 00 (...03) in the SNP.
  void compute_stats() {
    size_t nbByte = nbChars();
    size_t BitsInLastByte = (nbInds()%4); // number of bits to read on last byte

    /* FIRST : filling up stats_ with N0s, N1s, N2s, and NAs with PLINK translation*/

    // restarting with blank stats_ :
    stats_[0] = stats_[1] = stats_[2] = stats_[3] = 0;

    for (size_t i = 0; i < nbByte; i++) {
        uint8_t d = data()[i];
        if ((i == nbByte-1) && BitsInLastByte > 0) {
          while (BitsInLastByte > 0) {
            unsigned int val = (d & 3);
            unsigned int val_plink = currentMode_[0][val];
            stats_[val_plink]++;
            BitsInLastByte--;
            d >>= 2; // 1 shift par loop
          }
          if (!mu_ || !sigma_) compute_mu_sigma(); //compute mu & sd ONLY IF == 0 
          return compute_mode();
        }
      //sinon, fini pile à la fin d'une byte je peux return
      else if (i == nbByte-1) {
        if (!mu_ || !sigma_) compute_mu_sigma();
        return compute_mode();
      }
      stats_[0] += N0[d];
      stats_[2] += N0[255-d]; // raw = 3
      stats_[3] += N1[d]; // raw = 1
      stats_[1] += N1[255-d]; // raw = 2
    }
  }

  
  // for scalar product :
  // TODO : think on what are the limits (error checking) => what calls in gaston ?
  double LD(const SNPvector & other) {
    //if (nbInds() != other.nbInds()) // CAUSES PROBLEMS WITH CONST
    // TODO : check that read using same mode ?
    double LD = 0;
    double gg[16];
    gg[1] = gg[4] = gg[5] = gg[6] = gg[7] = gg[9] = gg[13] = 0;
    
    gg[0] = (-mu_)*(-(other.mu_));
    gg[2] = (-mu_)*(1.-(other.mu_));
    gg[3] = (-mu_)*(2.-(other.mu_));

    gg[8] = (1.-mu_)*(-(other.mu_));
    gg[10] = (1.-mu_)*(1.-(other.mu_));
    gg[11] = (1.-mu_)*(2.-(other.mu_));

    gg[12] = (2.-mu_)*(-(other.mu_));
    gg[14] = (2.-mu_)*(1.-(other.mu_));
    gg[15]= (2.-mu_)*(2.-(other.mu_));

    double v = sigma_ * other.sigma_;

#pragma omp parallel for reduction(+:LD)
    for(size_t i = 0; i < nbChars(); i++) {
      uint8_t g1 = data()[i]; //je récup les ièmes char
      const uint8_t g2_const = other.data()[i];
      uint8_t g2 = g2_const;
      for(int ss = 0; ss < 4; ss++) { // que je vais lire 2bits par 2bits
        // NOT using mode here cos gg modified to fit
        int g1val = g1&3;
        int g2val = g2&3;
        LD += gg[ (g1val*4) + g2val ];
        g1 >>= 2;
        g2 >>= 2;
      }
    }
    double r = LD/(v*(nbInds()-1)); //nbInds should be the same for both
    return r * r;
  }

  //function returning a 16 int vector holding frequency distribution
  std::vector<int> contingency(const SNPvector & other) {
    std::vector<int> table(16);
// #pragma omp declare reduction(vec_incr : std::vector<int> : \
//   std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
//   initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
// #pragma omp parallel for reduction(vec_incr : table)
    for(size_t i = 0; i < nbChars(); i++) {
      uint8_t snp1 = data()[i]; //je récup les ièmes char
      const uint8_t snp2_const = other.data()[i];
      uint8_t snp2 = snp2_const;
      for(int ss = 0; ss < 4; ss++) { // que je vais lire 2bits par 2bits
        // Hardcoded to only use Plink here, without touching mode
        int ind1 = currentMode_[0][snp1&3];
        int ind2 = currentMode_[0][snp2&3];
        table[ind1*4 + ind2]++;
        //std::cout << "ind1 = " << ind1 << " ind2 = " << ind2 << ", " << "Table[" << ind1*4 + ind2 << "] = " << table[ind1*4 + ind2] << "\n";
        ind1 >>= 2;
        ind2 >>= 2;
      }
    }
    return table;
  }
  


  class Iterator {
    private:
    int currentChar; // ii dans le code RV, correspond au byte sur lequel je suis
    int current2bits; // ss dans le code RV, correspond au bit (0...3) dans le byte
    SNPvector & iterated; //link to mother class

    public:
    // le constructeur vide démarre à 0, et récupère l'instance qui l'appelle
    Iterator(SNPvector & it) : currentChar(0), current2bits(0), iterated(it) {}

    Iterator(size_t ind, SNPvector & it) : currentChar(ind / 4), current2bits(ind % 4), iterated(it) {
      if (currentChar > iterated.nbChars()) throw std::out_of_range("End Iterator is too far !");
    }

    // operateur * const : renvoie la valeur 2bits par 2bits
    double operator*() {
      uint8_t byte = iterated.data()[currentChar];
      byte >>= (2 * current2bits); // pour avoir les 2bits de poids faible
      unsigned int val = (byte & 3); // mtn on les isole
      return iterated.mode(val);
    }

    // Operateur ++ : passe à la valeur suivante, PRE-INCRÉMENTATION 
    Iterator & operator++() {
      current2bits++;
      if (current2bits > 3) 
      {
        current2bits = 0;
        currentChar++;
      }
      return *this;
    }

    // Opérateur de comparaison entre deux itérateurs
    bool operator!=(const Iterator & other) const {
      return (currentChar != other.currentChar) || (current2bits != other.current2bits);
    }
  };

  // begin() : renvoie un itérateur
  Iterator begin() {
    return Iterator(*this);
  }

  /* end() : renvoie un itérateur qui "pointe vers la fin"
  aka après le dernier individu*/
  Iterator end() {
    return Iterator(nbInds(),*this);
  }
};

#endif