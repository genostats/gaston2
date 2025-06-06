#include <Rcpp.h>

#include <iostream>

#include "LD.h"
#include "SNPmatrix.h"
#include "mmatrix_methods.h"

// [[Rcpp::export]]
void LD_square_bigmemory(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  if (usefloat) {
    MMatrix<float> LD(path, i2 - i1 + 1, i2 - i1 + 1);
    LD_matrix<LDalgorithm::moments, float>(*pM, i1, i2, LD);
    LD.create_descriptor_file();
  } else {
    MMatrix<double> LD(path, i2 - i1 + 1, i2 - i1 + 1);
    LD_matrix<LDalgorithm::moments, double>(*pM, i1, i2, LD);
    LD.create_descriptor_file();
  }
}

// [[Rcpp::export]]
void LD_chunk_bigmemory(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2, size_t j1, size_t j2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  if (usefloat) {
    MMatrix<float> LD(path, i2 - i1 + 1, j2 - j1 + 1);
    LD_chunk<LDalgorithm::moments, float>(*pM, i1, i2, j1, j2, LD);
    LD.create_descriptor_file();
  } else {
    MMatrix<double> LD(path, i2 - i1 + 1, j2 - j1 + 1);
    LD_chunk<LDalgorithm::moments, double>(*pM, i1, i2, j1, j2, LD);
    LD.create_descriptor_file();
  }
}

// [[Rcpp::export]]
void LD_square_EM_bigmemory(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  if (usefloat) {
    MMatrix<float> LD(path, i2 - i1 + 1, i2 - i1 + 1);
    LD_matrix<LDalgorithm::EM, float>(*pM, i1, i2, LD);
    LD.create_descriptor_file();
  } else {
    MMatrix<double> LD(path, i2 - i1 + 1, i2 - i1 + 1);
    LD_matrix<LDalgorithm::EM, double>(*pM, i1, i2, LD);
    LD.create_descriptor_file();
  }
}

// [[Rcpp::export]]
void LD_chunk_EM_bigmemory(Rcpp::XPtr<SNPmatrix> pM, size_t i1, size_t i2, size_t j1, size_t j2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  if (usefloat) {
    MMatrix<float> LD(path, i2 - i1 + 1, j2 - j1 + 1);
    LD_chunk<LDalgorithm::EM, float>(*pM, i1, i2, j1, j2, LD);
    LD.create_descriptor_file();
  } else {
    MMatrix<double> LD(path, i2 - i1 + 1, j2 - j1 + 1);
    LD_chunk<LDalgorithm::EM, double>(*pM, i1, i2, j1, j2, LD);
    LD.create_descriptor_file();
  }
}