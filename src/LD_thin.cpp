#include "SNPmatrix.h"
#include "LD_thin.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
void LD_thin_(Rcpp::XPtr<SNPmatrix<>> pM, double threshold, int max_dist_bp, double max_dist_cM, Rcpp::LogicalVector which_keep)  {
  pM->computeSNPStats();
  LD_thin(*pM, threshold, max_dist_bp, max_dist_cM, which_keep, [](size_t i, size_t k) {return true;});
} 

