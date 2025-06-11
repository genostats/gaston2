#include <Rcpp.h>

#include "SNPmatrix.h"
#include "SNPdosageMemory.h"  //should include also SNPdosage

/* Function that eats a SNPmatrix<SNPdosage>
and copies it back it to SNPmatrix<SNPdosageMemory> specifically
with all the data and metadata from the original matrix.
If the original Matrix was on disk, this means putting all data in memory !
*/
void ToDosagematrixMemory(const SNPmatrix<SNPdosage> &other, SNPmatrix<SNPdosageMemory> &newMat) {

  const std::vector<std::shared_ptr<SNPdosage>> otherSNPs = other.getSNPs();
  for (const auto &snp : otherSNPs) {
    newMat.push_back(std::make_shared<SNPdosageMemory>(snp));
  }
  // extract stats now and set stats_set_ to true
  newMat.setIndStats(other.getIndStats());
  newMat.setSnpStats(other.getSNPStats());
}

SNPmatrix<SNPdosageMemory> ToDosagematrixMemory(SNPmatrix<SNPdosage> other, std::string newfile_name) {
    SNPmatrix<SNPdosageMemory> M;
    ToDosagematrixMemory(other, M);
    return M;
  }

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosageMemory>> ToDosagematrixMemory_(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM) {
  Rcpp::XPtr<SNPmatrix<SNPdosageMemory>> pnewMat(new SNPmatrix<SNPdosageMemory>);
  ToDosagematrixMemory(*pM,*pnewMat);
  return pnewMat;
}