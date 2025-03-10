#include <cstdint>
#include <cstddef>

//TODO cleanup if not useful
#include <iostream> // only for printing in debug

#ifndef _SNPvector_
#define _SNPvector_

  //TODO : find a way to change default value
  static uint8_t mu = 0;
  static uint8_t sigma = 1;

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
  
  enum Mode { 
    NUMERIC = 0, // (g = {0, 1, 2 and 3 = NA})
    CENTERED = 1, // (g -= mu)
    STANDARDIZED = 2, // (g = (g-mu)/sd)
    PLINK = 3// .bed file : (g = {0, 1, 3 and 2 = NA})
  };
  //TODO : check if done correctly 
  virtual Mode mode() = 0;
  //TODO : add N0, N1, N2, N3 = nbre d’élts de valeurs 0 1 2 et 3

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
          S += val;
          cptr++;
      }
    }
    return S;
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
      //if (val != 3) 
      //std::cout << "This is the " << currentChar << "th byte on the " << current2bits << " bit value :" << val << "\n";
      //std::cout << "This is the iterated.mode" << iterated.mode() << "\n";
      switch (iterated.mode()) {
            case NUMERIC:
                return val;
            case CENTERED:
                return (val - mu);
            case STANDARDIZED:
                return (val - mu) / sigma; // TODO : check if sigma CAN divide 
            case PLINK:
                if (val == 3) return 2;
                if (val == 2) return 3;
                return val;
            default:
                return val;
        }
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