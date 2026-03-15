#include "SNPmatrix.h"
#include "LD.h"
#include <iostream>
#include <Rcpp.h>

// TODO better handling of computation of stats ?!
// (if stats are already computed, these functions do nothing, but there's still
//  a function call and a test...)

// moments estimates. All in double

// [[Rcpp::export]]
double LD_pair_moments(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, bool r_scale) {
  pM->getSNP(i1)->compute_stats();
  pM->getSNP(i2)->compute_stats();
  auto f = LD_pair_f<LDalgorithm::moments>{};
  return f(*pM, i1, i2, r_scale);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix LD_square_moments(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, bool r_scale) {
  pM->computeSNPStats(i1, i2);
  Rcpp::NumericMatrix LD(i2 - i1 + 1, i2 - i1 + 1);
  LD_matrix<LDalgorithm::moments>(*pM, i1, i2, LD, r_scale);
  return LD;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix LD_chunk_moments(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, size_t j1, size_t j2, bool r_scale) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::NumericMatrix LD(i2 - i1 + 1, j2 - j1 + 1);
  LD_chunk<LDalgorithm::moments>(*pM, i1, i2, j1, j2, LD, r_scale);
  return LD;
}

// ----------------------------------------------------------------
// EM estimates. All in double

// [[Rcpp::export]]
double LD_pair_EM(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, bool r_scale) {
  pM->getSNP(i1)->compute_stats();
  pM->getSNP(i2)->compute_stats();
  auto f = LD_pair_f<LDalgorithm::EM>{};
  return f(*pM, i1, i2, r_scale);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix LD_square_EM(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, bool r_scale) {
  pM->computeSNPStats(i1, i2);
  Rcpp::NumericMatrix LD(i2 - i1 + 1, i2 - i1 + 1);
  LD_matrix<LDalgorithm::EM>(*pM, i1, i2, LD, r_scale);
  return LD;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix LD_chunk_EM(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, size_t j1, size_t j2, bool r_scale) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::NumericMatrix LD(i2 - i1 + 1, j2 - j1 + 1);
  LD_chunk<LDalgorithm::EM>(*pM, i1, i2, j1, j2, LD, r_scale);
  return LD;
}
