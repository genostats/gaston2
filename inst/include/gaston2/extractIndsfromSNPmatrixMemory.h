#include "SNPmatrix.h"
#include "SNPvectorMemory.h"
#include <iostream>
#include <memory> // shared_ptr
#include "mode.h" // for mode manipulation, tmp
#include <cstring>
#include <Rcpp.h>

#ifndef _EXTRACTINDMATRIX_MEM_
#define _EXTRACTINDMATRIX_MEM_


// filling up a new SNPmatrix with only the individuals at index specified in "keep"
template <typename SNPvectorClass, typename intVec>
void extractIndsfromSNPmatrixMemory(const SNPmatrix<SNPvectorClass> &other, const intVec &keep, SNPmatrix<SNPvectorClass> & newMat) {

  const std::vector<std::shared_ptr<SNPvector>> otherSNPs = other.getSNPs();

  for (const auto &snp : otherSNPs){
    newMat.push_back(std::make_shared<SNPvectorMemory>(snp, keep), false);
  }
  //extract stats now and set stats_set_ to true
  DataStruct original_dt = other.getIndStats();
  newMat.setIndStats(DataStruct(original_dt, keep));

  // keeping all SNPStats from file, but need to update the N0...
  newMat.setSnpStats(other.getSNPStats());

  newMat.exportSNPStats(true); // to force recompute and not just copy og N0...
  newMat.setMode(other.getMode());
}

template <typename SNPvectorClass, typename intVec>
SNPmatrix<SNPvectorClass> extractIndsfromSNPmatrixMemory(const SNPmatrix<SNPvectorClass> &other, const intVec &keep) {
  SNPmatrix<SNPvectorClass> M;
  extractIndsfromSNPmatrixMemory(other, keep, M);
  return M;
}

#endif 
