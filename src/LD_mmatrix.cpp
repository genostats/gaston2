#include <Rcpp.h>

#include "LD.h"
#include "SNPmatrix.h"
#include "houba/MMatrix.h"

// Fonctions à destination de la classe mmatrix en R

// [[Rcpp::export]]
SEXP LD_square_moments_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, std::string path, bool r_scale, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, i2 - i1 + 1);
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>> LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::moments, float>(*pM, i1, i2, *LD_ptr, r_scale);
    LD.slot("ptr") = LD_ptr;
    LD.slot("datatype") = "float";
  } else {
    Rcpp::XPtr<houba::MMatrix<double>> LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::moments, double>(*pM, i1, i2, *LD_ptr, r_scale);
    LD.slot("ptr") = LD_ptr;
    LD.slot("datatype") = "double";
  }
  return LD;
}

// [[Rcpp::export]]
SEXP LD_chunk_moments_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, size_t j1, size_t j2, std::string path, 
                              bool r_scale, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, j2 - j1 + 1);
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>> LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::moments, float>(*pM, i1, i2, j1, j2, *LD_ptr, r_scale);
    LD.slot("datatype") = "float";
    LD.slot("ptr") = LD_ptr;
  } else {
    Rcpp::XPtr<houba::MMatrix<double>> LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::moments, double>(*pM, i1, i2, j1, j2, *LD_ptr, r_scale);
    LD.slot("datatype") = "double";
    LD.slot("ptr") = LD_ptr;
  }
  return LD;
}

// [[Rcpp::export]]
SEXP LD_square_EM_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, std::string path, bool r_scale, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, i2 - i1 + 1);
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>> LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::EM, float>(*pM, i1, i2, *LD_ptr, r_scale);
    LD.slot("datatype") = "float";
    LD.slot("ptr") = LD_ptr;
  } else {
    Rcpp::XPtr<houba::MMatrix<double>> LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::EM, double>(*pM, i1, i2, *LD_ptr, r_scale);
    LD.slot("datatype") = "double";
    LD.slot("ptr") = LD_ptr;
  }
  return LD;
}

// [[Rcpp::export]]
SEXP LD_chunk_EM_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, size_t j1, size_t j2, std::string path,
                         bool r_scale, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, j2 - j1 + 1);
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>> LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::EM, float>(*pM, i1, i2, j1, j2, *LD_ptr, r_scale);
    LD.slot("datatype") = "float";
    LD.slot("ptr") = LD_ptr;
  } else {
    Rcpp::XPtr<houba::MMatrix<double>> LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::EM, double>(*pM, i1, i2, j1, j2, *LD_ptr, r_scale);
    LD.slot("datatype") = "double";
    LD.slot("ptr") = LD_ptr;
  }
  return LD;
}
