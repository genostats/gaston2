#include "SNPmatrix.h"
#include "SNPvectorMemory.h"
#include <iostream>
#include <memory> // shared_ptr
#include <cstring>
#include <Rcpp.h>

#ifndef _BINDINDMATRIX_MEM_
#define _BINDINDMATRIX_MEM_


// filling up a new SNPmatrix with every individuals from the first and second matrix
void bindIndstoSNPmatrixMemory(const SNPmatrix<SNPvector> &first, const SNPmatrix<SNPvector> &second, SNPmatrix<SNPvectorMemory> &newMat) {

  if (first.size() != second.size())
  throw std::logic_error("You cannot merge 2 SNPmatrix with mismatched number of SNPs (maybe later could add NAs)");

  const std::vector<std::shared_ptr<SNPvector>> firstSNPs = first.getSNPs();
  const std::vector<std::shared_ptr<SNPvector>> secondSNPs = second.getSNPs();

  // maybe could also do a constructor or a function on the SNPmatrix level ?
  for (size_t i = 0; i < firstSNPs.size(); i++) {
    newMat.push_back(std::make_shared<SNPvectorMemory>(firstSNPs[i], secondSNPs[i]));
  }

  // Individuals stats are still good (N0, N1, N2, NAs also)
  // so I can just fuse them
  newMat.setIndStats(DataStruct(first.getIndStats(), second.getIndStats()));
  // SNP stats need to be recomputed
  newMat.setSnpStats(first.getSNPStats());

}

SNPmatrix<SNPvectorMemory> bindIndstoSNPmatrixMemory(const SNPmatrix<SNPvector> &first, const SNPmatrix<SNPvector> &second) {
  SNPmatrix<SNPvectorMemory> M;
  bindIndstoSNPmatrixMemory(first, second, M);
  return M;
}

#endif