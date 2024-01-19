#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector>
#include <memory>

#ifndef _snpmatrix_
#define _snpmatrix_

/**
 * @brief A class keeping in a vector pointers to SNPvectors
 * (either SNPvectorDisk or SNPvectorMemory)
 * 
 * and not a very good idea to leave SNPs public
 * --> will change
 */
class SNPmatrix {
  public:

  void push_back(std::shared_ptr<SNPvector> v) {
    SNPs.push_back(v);
  }

  int size() { return SNPs.size(); }

  //temporary func to test d°
  void deleteSNP() {
    SNPs.pop_back();
  }

  /**
   * @brief A member function testing whether the SNP at index 
   * is a SNPVectorDisk or a SNPVectorMemory
   * 
   * @param index 
   * @return true 
   * @return false 
   */
  bool onDisk(size_t index) {
    if (std::shared_ptr<SNPVectorDisk> test = std::dynamic_pointer_cast<SNPVectorDisk>(SNPs[index])) return true;
    return false;
  }

  std::vector<std::shared_ptr<SNPvector>> SNPs;
};

#endif