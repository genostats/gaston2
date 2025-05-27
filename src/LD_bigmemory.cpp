#include "SNPmatrix.h"
#include "LD.h"
#include <iostream>
#include <Rcpp.h>
#include "mmatrix_methods.h"

// [[Rcpp::export]]
void LD_square_bigmemory(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2, std::string path) {
  pM->computeSNPStats(i1, i2);
  MMatrix<float> LD(path, i2 - i1 + 1, i2 - i1 + 1);
  LD_matrix<LDalgorithm::moments, float>(*pM, i1, i2, LD);
  LD.create_descriptor_file();
}

/*
// [[Rcpp::export]]
Rcpp::NumericMatrix LD_chunk(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2, size_t j1, size_t j2) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::NumericMatrix LD(i2 - i1 + 1, j2 - j1 + 1);
  LD_chunk<LDalgorithm::moments>(*pM, i1, i2, j1, j2, LD);
  return LD;
}
*/
