#include "SNPmatrix.h"
#include "SNPdosage.h"
#include "mode.h"
#include <Rcpp.h>


// [[Rcpp::export]]
void computeIndStats(Rcpp::XPtr<SNPmatrix<>> pM, bool force = false) {
  pM->compute_indStats(force);
}

// a priori inutile ? la première fonction doit marcher quelque soit le cas
// [[Rcpp::export]]
std::string computeIndStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool force = false) {
  pM->compute_indStats(force);
}
