#include <cstdint>
#include <cstddef>

//TODO cleanup if not usefull
#include <iostream> // only for printing in debug

#ifndef _SNPvector_
#define _SNPvector_
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

  // // an idiotic function just for testing purposes
  // // will allow us to compare performances
  // unsigned int sum() {
  //   unsigned int S = 0;
  //   uint8_t * d = data();
  //   size_t nc = nbChars();
  //   for(size_t i = 0; i < nc; i++) {
  //     S += (unsigned int) d[i];
  //   }
  //   return S;
  // }

  unsigned int sum() {
    unsigned int S = 0;
    uint8_t * d = data();
    size_t nc = nbChars();
    
    for(size_t i = 0; i < nc; i++) {
      uint8_t byte = d[i];
      for (unsigned int bits = 0; bits < 4; bits++) {
        byte >>= (2 * bits); // pour avoir les 2bits de poids faible
        S += (byte & 3);
      }
    }
    
    return S;
  }

  // TODO : check if class Iterator is inherited to SNPvectorDisk and SNPvectorMemory
  class Iterator {
    private:
    int currentChar; // ii dans le code RV, correspond au byte sur lequel je suis
    int current2bits; // ss dans le code RV, correspond au bit (0...3) dans le byte

    uint8_t * iterated;

    public:
    // le constructeur vide démarre à 0, et récupère l'instance qui l'appelle
    Iterator(uint8_t * it) : currentChar(0), current2bits(0), iterated(it) {}

    // TODO : bound check to add
    Iterator(size_t ind, uint8_t * it) : currentChar(ind / 4), current2bits(ind % 4), iterated(it) {}

    // operateur * const : renvoie la valeur 2bits par 2bits
    unsigned int operator*() {
      uint8_t byte = iterated[currentChar];
      byte >>= (2 * current2bits); // pour avoir les 2bits de poids faible
      unsigned int val = (byte & 3); // mtn on les isole
      // if (val != 3) std::cout << "This is the " << currentChar << "th byte on the " << current2bits << " bit value :" << val << "\n";
      return val;
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
    return Iterator(data());
  }

  /* end() : renvoie un itérateur qui "pointe vers la fin"
  aka  vers le dernier individu*/
  Iterator end() {
    return Iterator(nbInds(), data());
  }
};

#endif