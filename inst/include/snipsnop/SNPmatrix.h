#include "SNPvector.h"
#include <vector>
#include <memory>

#ifndef _snpmatrix_
#define _snpmatrix_

// not much yet here 
// and not a very good idea to leave SNPs public
// --> will change
class SNPmatrix {
  public:
  std::vector<std::shared_ptr<SNPvector>> SNPs;
  void push_back(std::shared_ptr<SNPvector> v) {
    SNPs.push_back(v);
  }
};


#endif

