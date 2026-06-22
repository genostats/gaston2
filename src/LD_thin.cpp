#include "SNPmatrix.h"
#include "LD_thin.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
void LD_thin_left(Rcpp::XPtr<SNPmatrix<>> pM, double threshold, int max_dist_bp, double max_dist_cM, Rcpp::LogicalVector which_keep)  {
  pM->computeSNPStats();
  LD_thin(*pM, threshold, max_dist_bp, max_dist_cM, which_keep, [](size_t i, size_t k) {return true;});
} 

// [[Rcpp::export]]
void LD_thin_right(Rcpp::XPtr<SNPmatrix<>> pM, double threshold, int max_dist_bp, double max_dist_cM, Rcpp::LogicalVector which_keep)  {
  pM->computeSNPStats();
  LD_thin(*pM, threshold, max_dist_bp, max_dist_cM, which_keep, [](size_t i, size_t k) {return false;});
}

// [[Rcpp::export]]
void LD_thin_random(Rcpp::XPtr<SNPmatrix<>> pM, double threshold, int max_dist_bp, double max_dist_cM, Rcpp::LogicalVector which_keep) {
  pM->computeSNPStats();
  LD_thin(*pM, threshold, max_dist_bp, max_dist_cM, which_keep, [](size_t i, size_t j) { return R::unif_rand() > 0.5; });
}

// [[Rcpp::export]]
void LD_thin_priority(Rcpp::XPtr<SNPmatrix<>> pM, double threshold, int max_dist_bp, double max_dist_cM, Rcpp::LogicalVector which_keep, 
                       Rcpp::NumericVector priority) {
  pM->computeSNPStats();
  LD_thin(*pM, threshold, max_dist_bp, max_dist_cM, which_keep, [&priority](size_t i, size_t j) { return priority[i] > priority[j]; });
}

