#include "SNPmatrix.h"
#include "SNPdosage.h"
#include <Rcpp.h>

// [[Rcpp::export]]
void computeSNPStats(Rcpp::XPtr<SNPmatrix<SNPvector>> pM) {
  pM->computeSNPStats();
}

// [[Rcpp::export]]
void computeSNPStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM) {
  pM->computeSNPStats();
}

// [[Rcpp::export]]
void exportSNPStats(Rcpp::XPtr<SNPmatrix<SNPvector>> pM, bool force = false) {
  pM->exportSNPStats(force);
}

// [[Rcpp::export]]
void exportSNPStatsDosage(Rcpp::XPtr<SNPmatrix<SNPdosage>> pM, bool force = false) {
  pM->exportSNPStats(force);
}
