#include <Rcpp.h>

#include "SNPmatrix.h"
#include "SNPvectorMemory.h"  //should include also SNPvector

/* Function that eats a SNPmatrix<SNPvector>
and copies it back it to SNPmatrix<SNPvectorMemory> specifically
with all the data and metadata from the original matrix.
If the original Matrix was on disk, this means putting all data in memory !
*/
void ToSNPmatrixMemory(const SNPmatrix<SNPvector> &other, SNPmatrix<SNPvectorMemory> &newMat) {

  const std::vector<std::shared_ptr<SNPvector>> otherSNPs = other.getSNPs();
  for (const auto &snp : otherSNPs) {
    newMat.push_back(std::make_shared<SNPvectorMemory>(snp));
  }
  // extract stats now and set stats_set_ to true
  newMat.setIndStats(other.getIndStats());
  newMat.setSnpStats(other.getSNPStats());
}

SNPmatrix<SNPvectorMemory> ToSNPmatrixMemory(SNPmatrix<> other, std::string newfile_name) {
    SNPmatrix<SNPvectorMemory> M;
    ToSNPmatrixMemory(other, M);
    return M;
  }

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPvectorMemory>> ToSNPmatrixMemory_(Rcpp::XPtr<SNPmatrix<SNPvector>> pM) {
  Rcpp::XPtr<SNPmatrix<SNPvectorMemory>> pnewMat(new SNPmatrix<SNPvectorMemory>);
  ToSNPmatrixMemory(*pM,*pnewMat);
  return pnewMat;
}