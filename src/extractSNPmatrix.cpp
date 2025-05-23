/**
 * @file extractSNPmatrix.cpp
 * @author Hervé Perdry
 * @brief Main file with functions to interact with SNPvector and SNPmatrix
 * exported by readBedFileMemory.h and called in test_readBedFileMemory.cpp
 * @date 2023-12-11
 * 
 */
#include "SNPmatrix.h"
#include "SNPvectorMemory.h"
#include <iostream>
#include <memory> // shared_ptr
#include "mio.hpp"
#include "SNPvectorDisk.h"
#include <cstring>
#include <Rcpp.h>
#include "extractSNPmatrix.h"

// TODO : think, do I want to return whole matrix ? ptr to matrix ?

//R exported function, keep is an Integer Vector by default, hope that's okay ?
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix> extractSNPmatrixMemory_(Rcpp::XPtr<SNPmatrix> other, Rcpp::IntegerVector &keep) {
  Rcpp::XPtr<SNPmatrix> pnewMat(new SNPmatrix);
  extractSNPmatrixMemory(*other, keep, *pnewMat);
  return pnewMat;
}

//R exported function, keep is an Integer Vector by default, hope that's okay ?
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix> extractSNPmatrixDisk_(Rcpp::XPtr<SNPmatrix> other, Rcpp::IntegerVector &keep, std::string path_str) {
  Rcpp::XPtr<SNPmatrix> pnewMat(new SNPmatrix);
  extractSNPmatrixDisk(*other, keep, path_str, *pnewMat);
  return pnewMat;
}