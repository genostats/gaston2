#include "SNPmatrix.h"
#include "SNPdosageMemory.h" // will include SNPdosage also
#include <iostream>
#include <memory> // shared_ptr
#include <cstring>
#include <Rcpp.h>

#ifndef _EXTRACTINDDOSMATRIX_MEM_
#define _EXTRACTINDDOSMATRIX_MEM_


// filling up a new SNPmatrix with only the individuals at index specified in "keep"
template <typename SNPvectorClass, typename intVec>
void extractIndsfromDosagematrixMemory(const SNPmatrix<SNPvectorClass> &other, const intVec &keep, SNPmatrix<SNPvectorClass> & newMat) {

  const std::vector<std::shared_ptr<SNPdosage>> otherSNPs = other.getSNPs();

  for (const auto &snp : otherSNPs){
    newMat.push_back(std::make_shared<SNPdosageMemory>(snp, keep));
  }
  // extract stats now and set stats_set_ to true
  DataStruct original_dt = other.getIndStats();
  newMat.setIndStats(DataStruct(original_dt, keep));
  // keeping all SNPStats
  newMat.setSnpStats(other.getSNPStats());

  newMat.computeSNPStats();
  newMat.setMode(other.getMode());
}

template <typename SNPvectorClass, typename intVec>
SNPmatrix<SNPvectorClass> extractIndsfromDosagematrixMemory(const SNPmatrix<SNPvectorClass> &other, const intVec &keep) {
  SNPmatrix<SNPvectorClass> M;
  extractIndsfromDosagematrixMemory(other, keep, M);
  return M;
}

#endif 
