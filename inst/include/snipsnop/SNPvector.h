#include <cstdint>
#include <cstddef>

//TODO cleanup if not useful
#include <iostream> // only for printing in debug

#ifndef _SNPvector_
#define _SNPvector_

  //TODO : find a way to change default value
  static uint8_t mu = 0;
  static uint8_t sigma = 1;

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

  // BY DEFAULT, PLINK format, with 01 = missing genotype
  constexpr uint8_t Defaultmode[4] = {0, 3, 1, 2};

/**
 * @brief An abstract class instanciated through SNPvectorMemory or SNPvectorDisk
 * 
 * Stores a bit vector, could be in memory or in a memory_mapped file
 * 
 */
class SNPvector {
  public:
  /**
   * @brief pure virtual function, 
   * returning a pointer to the front of the vector
   * 
   * @return uint8_t* 
   */
  virtual uint8_t * data() = 0;
  // nombre d'individus
  virtual size_t nbInds() = 0;

  virtual const uint8_t* mode() = 0;
  virtual uint8_t mode(unsigned int n) = 0;
  
  // TODO : maybe with env var ?
  //{0, 3, 1, 2};

  // TODO : A REFAIRE 
  // enum Mode { 
  //   NUMERIC = 0, // (g = {0, 1, 2 and 3 = NA})
  //   CENTERED = 1, // (g -= mu)
  //   STANDARDIZED = 2, // (g = (g-mu)/sd)
  //   PLINK = 3// .bed file : (g = {0, 1, 3 and 2 = NA})
  // };
  //TODO : check if done correctly 
  //virtual Mode mode() = 0;
  /* Containing N0, N1, N2, N3 on the whole SNP
  populated by compute_stats() */
  unsigned int stats[4] = {0, 0, 0, 0};

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

  // dummy function summing the n first individuals from the SNP
  unsigned int sum(int n = -1) {
    unsigned int S = 0;
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

  // Method filling up stats[] w/ the nb of ind = 00 (...03) in the SNP.
  // TODO : REVOIR POUR CLEAN CETTE METHODE TRES MOCHE ààààà
  unsigned int * compute_stats() {
    size_t nbByte = nbChars();
    size_t BitsInLastByte = (nbInds()%4); // number of bits to read on last byte
    for (size_t i = 0; i < nbByte; i++) {
        
        uint8_t d = data()[i];

        if ((i == nbByte-1) && BitsInLastByte > 0) {
          while (BitsInLastByte > 0) {
          
            //std::cout << "On last byte, " << BitsInLastByte << " left\n";
            unsigned int val = (d & 3);
            //std::cout << "This is the " << i << "th byte on the " << 3 - BitsInLastByte << " bit value :" << val << "\n";
            //std::cout << "Stats array before adding val :" << stats[0] << "," << stats[1]<< "," << stats[2]<< "," << stats[3] << "\n";
            stats[val]++;
            BitsInLastByte--;
            d >>= 2; // 1 shift par loop
          }
          return stats;
        }
      //sinon, fini pile à la fin d'une byte je peux return
      else if (i == nbByte-1) return stats;
      stats[0] += N0[d];
      stats[3] += N0[255-d];
      stats[1] += N1[d];
      stats[2] += N1[255-d];
    }
    // maybe add a translation following diff modes ? 
    return stats;
  }

  class Iterator {
    private:
    int currentChar; // ii dans le code RV, correspond au byte sur lequel je suis
    int current2bits; // ss dans le code RV, correspond au bit (0...3) dans le byte
    SNPvector & iterated; //link to mother class

    public:
    // le constructeur vide démarre à 0, et récupère l'instance qui l'appelle
    Iterator(SNPvector & it) : currentChar(0), current2bits(0), iterated(it) {}

    // TODO : bound check to add
    Iterator(size_t ind, SNPvector & it) : currentChar(ind / 4), current2bits(ind % 4), iterated(it) {
    //std::cout << "Creating an End iterator at " << currentChar << "th byte on the " << current2bits << " bit \n";
    }

    // operateur * const : renvoie la valeur 2bits par 2bits
    unsigned int operator*() {
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
      // TODO : modify with < maybe ?
      return (currentChar != other.currentChar) || (current2bits != other.current2bits);
    }
  };

  //begin et end dans la classe mère.

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