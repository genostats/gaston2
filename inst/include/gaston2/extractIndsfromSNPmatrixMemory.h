#include "SNPmatrix.h"
#include "SNPvectorMemory.h"
#include <iostream>
#include <memory> // shared_ptr
#include <cstring>
#include <Rcpp.h>

#ifndef _EXTRACTINDMATRIX_MEM_
#define _EXTRACTINDMATRIX_MEM_


// filling up a new SNPmatrix with only the individuals at index specified in "keep"
template <typename intVec>
void extractIndsfromSNPmatrixMemory(const SNPmatrix &other, const intVec &keep, SNPmatrix & newMat) {

  const std::vector<std::shared_ptr<SNPvector>> otherSNPs = other.getSNPs();

  for (const auto &snp : otherSNPs){
    newMat.push_back(std::make_shared<SNPvectorMemory>(snp, keep));
  }
  //extract stats now and set stats_set_ to true
  DataStruct original_dt = other.getIndStats();
  newMat.setIndStats(DataStruct(original_dt, keep));
}

template <typename intVec>
SNPmatrix extractIndsfromSNPmatrixMemory(const SNPmatrix &other, const intVec &keep) {
  SNPmatrix M;
  extractIndsfromSNPmatrixMemory(other, keep, M);
  return M;
}

#endif 
