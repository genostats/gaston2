#include <cstdint>
#include <cstddef>

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

  // an idiotic function just for testing purposes
  // will allow us to compare performances
  unsigned int sum() {
    unsigned int S = 0;
    uint8_t * d = data();
    size_t nc = nbChars();
    for(size_t i = 0; i < nc; i++) {
      S += (unsigned int) d[i];
    }
    return S;
  }

  // TODO : check if class Iterator is inherited to SNPvectorDisk and SNPvectorMemory
  class Iterator {
    private:
    /* pour l'instant je fais de le choix de stocker dans currentValue 
    le compteur d'indice, transmis ensuite à data pour obtenir la valeur à cet indice*/ 
    int currentValue;

    SNPvector & iterated;

    public:
    // le constructeur vide démarre à 0, et récupère l'instance qui l'appelle
    Iterator(SNPvector & it) : currentValue(0), iterated(it) {
      //std::cout << "Creating an Iterator \n";
      }

    // TODO : bound check to add
    Iterator(size_t n, SNPvector & it) : currentValue(n), iterated(it) {}

    // operateur * const : renvoie la valeur DANS LE SNP
    // PB, Iterator has no reference to its mother class ????
    unsigned int operator*() {
      uint8_t * d = iterated.data();
      return (unsigned int) d[currentValue];
    }

    // Operateur ++ : passe à la valeur suivante, PRE-INCRÉMENTATION 
    Iterator & operator++() {
      currentValue++;
      return *this;
    }

    // Opérateur de comparaison entre deux itérateurs
    bool operator!=(const Iterator & other) const {
      return currentValue < other.currentValue;
    }
  };

  //begin et end dans la classe mère.

  // begin() : renvoie un itérateur
  Iterator begin() {
    return Iterator(*this);
  }

  /* end() : renvoie un itérateur qui "pointe vers la fin"
  aka  vers le dernier individu TODO : add -1 ou pas ?*/
  Iterator end() {
    return Iterator( nbChars(), *this);
  }

  // at the end so it can access Iterator's ctor
  // TODO : try without, don't see how can be usefull
  // just a dangling reference to something never called ?
  // private :
  // Iterator it { Iterator(*this)};
  
};

#endif