#include <cstdint>
#include <cstddef>
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
    size_t n = nbInds();
    return n/4 + ((n%4 == 0u)?0:1);
  }

  // an idiotic function just for testing purposes
  unsigned int sum() {
    unsigned int S = 0;
    uint8_t * d = data();
    size_t nc = nbChars();
    for(size_t i = 0; i < nc; i++) {
      S += (unsigned int) d[i];
    }
    return S;
  }
};

#endif