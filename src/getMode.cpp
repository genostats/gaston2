#include "SNPmatrix.h"
#include "SNPdosage.h"
#include "mode.h"
#include <Rcpp.h>


// [[Rcpp::export]]
std::string getMode(Rcpp::XPtr<SNPmatrix<>> pM) {
  std::cout << "Getting mode for snp.matrix object\n";
  return modeToString(pM->mode());
}


// [[Rcpp::export]]
std::string getModeDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM) {
  std::cout << "Getting mode for dose.matrix object\n";
  return modeToString(pM->mode());
}
