#include "SNPmatrix.h"
#include "LD.h"
#include <iostream>
#include <Rcpp.h>

// !! these functions send back D 

// [[Rcpp::export]]
double LD_pair_EM(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2) {
  pM->getSNP(i1)->compute_stats();
  pM->getSNP(i2)->compute_stats();
  auto f = LD_pair_f<LDalgorithm::EM>{};
  return f(*pM, i1, i2);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix LD_square_EM(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2) {
  pM->computeSNPStats(i1, i2);
  Rcpp::NumericMatrix LD(i2 - i1 + 1, i2 - i1 + 1);
  LD_matrix<LDalgorithm::EM>(*pM, i1, i2, LD);
  return LD;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix LD_chunk_EM(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2, size_t j1, size_t j2) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::NumericMatrix LD(i2 - i1 + 1, j2 - j1 + 1);
  LD_chunk<LDalgorithm::EM>(*pM, i1, i2, j1, j2, LD);
  return LD;
}


