/**
 * @file readBedFileMemory.cpp
 * @author Hervé Perdry
 * @brief Main file with functions to interact with SNPvector and SNPmatrix
 * exported by readBedFileMemory.h and called in test_readBedFileMemory.cpp
 * @date 2023-12-11
 * 
 */
#include "SNPmatrix.h"
#include "SNPvectorMemory.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <memory>
#include "mio.hpp"
#include "SNPvectorDisk.h"
#include <cstring>
#include <Rcpp.h>

/**
 * @brief Reading a bed file, storing SNPs in a SNPmatrix returned
 * 
 * @param filename The name of the file to be opened
 * @param n_ind The number of individuals/samples
 * @param n_snp The number of SNP to read from the file and to load into the Matrix.
 * 
 * @return SNPmatrix, stocking shared_ptrs to SNPvectorMemory 
 */

void readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp, SNPmatrix & M) {
  std::ifstream file(filename, std::ifstream::binary);
  if(!file.is_open()) {
    throw std::runtime_error("Cannot open file " + filename);
  }

  // check magic number
  char magic[3];
  file.read(magic, 3);
  if(magic[0] != 108 || magic[1] != 27) {
    throw std::runtime_error("Not a bed file");
  }
  if(magic[2] != 1) {
    throw std::runtime_error("Not a bed file in SNP major mode");
  } 

  for(size_t i = 0; i < n_snp; i++) {
    //makes a shared_ptr on a vector of snips 
    std::shared_ptr<SNPvectorMemory> snpVec(new SNPvectorMemory(n_ind));
    size_t n = snpVec->nbChars();
    uint8_t * data = snpVec->data();
    file.read(reinterpret_cast<char *>(data), n);
    M.push_back(snpVec);
  }
  
  file.close();
}

SNPmatrix readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp) {
  SNPmatrix M;
  readBedFileMemory(filename, n_ind, n_snp, M);
  return M;
}

// R exported function
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix> readBedFileMemory_(std::string filename, size_t n_ind, size_t n_snp) {
  Rcpp::XPtr<SNPmatrix> pM(new SNPmatrix);
  readBedFileMemory(filename, n_ind, n_snp, *pM);
  return pM;
}