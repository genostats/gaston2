#include <Rcpp.h>

#include "LD.h"
#include "SNPmatrix.h"
#include "houba/MMatrix.h"

// Fonctions à destination de la classe mmatrix en R, 
// Je ne sais pas à quel point elles seraient utiles / utilisées
// dans le code en C++
// donc j'ai décidé de ne pas me préoccuper de l'interface pour l'instant
// d'autant plus qu'elles ont pour intérêt de garder le mapping 
// fait pas MMatrix au lieu d'en réouvrir un, sinon juste des appels 
// à LD_chunk et LD_matrix

// sinon interfaçage possible et relativement simple avec 2 helpers 
// (une pour C++ une pour R, créant la MMatrix avec un Rcpp::XPtr ou un *) 
// et la fonction principale qui prends une MMatrix par référence et la modifie en place

// [[Rcpp::export]]
SEXP LD_square_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, i2 - i1 + 1);;
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>>  LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::moments, float>(*pM, i1, i2, *LD_ptr);
    LD.slot("ptr") = LD_ptr;
    LD.slot("datatype") = "float";
  } else {
    Rcpp::XPtr<houba::MMatrix<double>> LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::moments, double>(*pM, i1, i2, *LD_ptr);
    LD.slot("ptr") = LD_ptr;
    LD.slot("datatype") = "double";
  }
  return LD;
}

// [[Rcpp::export]]
SEXP LD_chunk_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, size_t j1, size_t j2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, j2 - j1 + 1);
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>>  LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::moments, float>(*pM, i1, i2, j1, j2, *LD_ptr);
    LD.slot("datatype") = "float";
    LD.slot("ptr") = LD_ptr;
  } else {
    Rcpp::XPtr<houba::MMatrix<double>>  LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::moments, double>(*pM, i1, i2, j1, j2, *LD_ptr);
    LD.slot("datatype") = "double";
    LD.slot("ptr") = LD_ptr;
  }
  return LD;
}

// [[Rcpp::export]]
SEXP LD_square_EM_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, i2 - i1 + 1);
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>>  LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::EM, float>(*pM, i1, i2, *LD_ptr);
    LD.slot("datatype") = "float";
    LD.slot("ptr") = LD_ptr;
  } else {
    Rcpp::XPtr<houba::MMatrix<double>>  LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, i2 - i1 + 1));
    LD_matrix<LDalgorithm::EM, double>(*pM, i1, i2, *LD_ptr);
    LD.slot("datatype") = "double";
    LD.slot("ptr") = LD_ptr;
  }
  return LD;
}

// [[Rcpp::export]]
SEXP LD_chunk_EM_mmatrix(Rcpp::XPtr<SNPmatrix<>> pM, size_t i1, size_t i2, size_t j1, size_t j2, std::string path, bool usefloat = true) {
  pM->computeSNPStats(i1, i2);
  pM->computeSNPStats(j1, j2);
  Rcpp::S4 LD("mmatrix");
  LD.slot("file") = path;
  LD.slot("dim") = Rcpp::IntegerVector::create(i2 - i1 + 1, j2 - j1 + 1);
  LD.slot("readonly") =  false;
  if (usefloat) {
    Rcpp::XPtr<houba::MMatrix<float>>  LD_ptr(new houba::MMatrix<float>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::EM, float>(*pM, i1, i2, j1, j2, *LD_ptr);
    LD.slot("datatype") = "float";
    LD.slot("ptr") = LD_ptr;
  } else {
    Rcpp::XPtr<houba::MMatrix<double>>  LD_ptr(new houba::MMatrix<double>(path, i2 - i1 + 1, j2 - j1 + 1));
    LD_chunk<LDalgorithm::EM, double>(*pM, i1, i2, j1, j2, *LD_ptr);
    LD.slot("datatype") = "double";
    LD.slot("ptr") = LD_ptr;
  }
  return LD;
}
