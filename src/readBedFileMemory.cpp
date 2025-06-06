/**
 * @file readBedFileMemory.cpp
 * @author Hervé Perdry
 * @brief Main file with functions to interact with SNPvector and SNPmatrix
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
 * @param bedfile The name of the bed file
 * @param bimfile The name of the bim file
 * @param bamfile The name of the fam file
 * @param M reference to an originally empty SNPmatrix
 */

void readBedFileMemory(std::string bedfile, std::string bimfile, std::string famfile, SNPmatrix<> & M) {
  // first read bim and fam file, to determine size
  M.readFamFile(famfile);
  M.readBimFile(bimfile);

  size_t n_snp = M.getSNPStats().nrow();
  size_t n_ind = M.getIndStats().nrow();

  // now ready to read bed file
  std::ifstream file(bedfile, std::ifstream::binary);
  if(!file.is_open()) {
    throw std::runtime_error("Cannot open file " + bedfile);
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

SNPmatrix<> readBedFileMemory(std::string bedfile, std::string bimfile, std::string famfile) {
  SNPmatrix<> M;
  readBedFileMemory(bedfile, bimfile, famfile, M);
  return M;
}

// R exported function
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<>> readBedFileMemory_(std::string bedfile, std::string bimfile, std::string famfile) {
  Rcpp::XPtr<SNPmatrix<>> pM(new SNPmatrix<>);
  readBedFileMemory(bedfile, bimfile, famfile, *pM);
  return pM;
}
