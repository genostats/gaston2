#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector>
#include <memory>

#ifndef _snpmatrix_
#define _snpmatrix_

// not much yet here 
// and not a very good idea to leave SNPs public
// --> will change
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

  bool onDisk(size_t index) { //doesn't work, will see later
    if (std::shared_ptr<SNPVectorDisk> test = std::dynamic_pointer_cast<SNPVectorDisk>(SNPs[index])) return true;
    return false;
  }

  std::vector<std::shared_ptr<SNPvector>> SNPs;
};

#endif