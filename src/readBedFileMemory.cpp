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
 * @param mode Mode index for reading the SNP, with NUMERIC = 0 (g = {0, 1, 2 and 3 = NA}), CENTERED = 1(g -= mu),
 * STANDARDIZED = 2 (g = (g-mu)/sd), PLINK = 3 .bed file : (g = {0, 1, 3 and 2 = NA})
 * 
 * @return SNPmatrix, stocking shared_ptrs to SNPvectorMemory 
 */
SNPmatrix readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp, Mode modeArray = PLINK) {
  std::ifstream file(filename, std::ifstream::binary);
  if(!file.is_open()) {
    throw std::runtime_error("Cannot open file");
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

  SNPmatrix M;
  for(size_t i = 0; i < n_snp; i++) {
    //makes a shared_ptr on a vector of snips 
    std::shared_ptr<SNPvectorMemory> snpVec(new SNPvectorMemory(n_ind, modeArray));
    size_t n = snpVec->nbChars();
    uint8_t * data = snpVec->data();
    file.read(reinterpret_cast<char *>(data), n);
    M.push_back(snpVec);
  }
  
  file.close();
  return M;
}

// R exported function
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix> readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp) {
std::ifstream file(filename, std::ifstream::binary);
  if(!file.is_open()) {
    throw std::runtime_error("Cannot open file");
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

  Rcpp::XPtr<SNPmatrix> pM(new SNPmatrix);
  for(size_t i = 0; i < n_snp; i++) {
    //makes a shared_ptr on a vector of snips 
    std::shared_ptr<SNPvectorMemory> snpVec(new SNPvectorMemory(n_ind));
    size_t n = snpVec->nbChars();
    uint8_t * data = snpVec->data();
    file.read(reinterpret_cast<char *>(data), n);
    pM->push_back(snpVec);
  }

  file.close();
  return pM;
}


