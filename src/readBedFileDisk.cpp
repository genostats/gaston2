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

void readBedFileDisk(std::string bedfile, std::string bimfile, std::string famfile, SNPmatrix<> & M) {
  // first read bim and fam file, to determine size
  M.readFamFile(famfile);
  M.readBimFile(bimfile);

  size_t n_snp = M.getSNPStats().nrow();
  size_t n_ind = M.getIndStats().nrow();

  // now ready to read bed file

  std::ifstream file_test(bedfile, std::ifstream::binary);
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  std::error_code error;
  mio::mmap_source file_ = mio::make_mmap_source(bedfile, 0, mio::map_entire_file, error);
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
  }
  std::shared_ptr<mio::mmap_source> file_ptr = std::make_shared<mio::mmap_source>(std::move(file_)); // don't know if good idea, creates a NEW sink
  // DON'T USE FILE_ FROM NOW ON COS IT'S NULL, WAS MOVED
  const char* data = reinterpret_cast<const char*>(file_ptr->data());// const necessary bcos read only 
  
  // check magic number
  char magic[3];
  for (int i = 0; i < 3; i++)
  {
    //std::cout << int(*data) << ' ';
    magic[i] = *(data + i);
  }
  //std::cout << '\n';
  if(magic[0] != 108 || magic[1] != 27) {
    throw std::runtime_error("Not a bed file");
  }
  if(magic[2] != 1) {
    throw std::runtime_error("Not a bed file in SNP major mode");
  }

  for(size_t i = 0; i < n_snp; i++) {
    std::shared_ptr<SNPvectorDisk<mio::access_mode::read>> snpVec(new SNPvectorDisk<mio::access_mode::read>(n_ind,file_ptr, i));
    //data should be taken by file_ptr

    size_t n = snpVec->nbChars(); // func inherited from SNPVec, gives back size used by SNP
    M.push_back(snpVec);
  }
}

/** @fn SNPmatrix readBedFileDisk(std::string path, size_t n_ind, size_t n_snp)
 * @brief Reading a bed file with memory mapping, storing SNPs in a SNPmatrix returned
 * 
 * A more detailled description here
 * 
 * @param path The path to the file to be opened. Can be relative or absolute
 * @param n_ind The number of individuals/samples, given by the .bim ?
 * @param n_snp The number of SNP to read from the file and to load into the Matrix.
 *   
 * @return a SNPmatrix, stocking shared_ptrs to SNPvectorDisk 
*/

SNPmatrix<> readBedFileDisk(std::string bedfile, std::string bimfile, std::string famfile) {
  SNPmatrix<> M;
  readBedFileDisk(bedfile, bimfile, famfile, M);
  return M;
}


// R exported function
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<>> readBedFileDisk_(std::string bedfile, std::string bimfile, std::string famfile) {
  Rcpp::XPtr<SNPmatrix<>> pM(new SNPmatrix<>);
  readBedFileDisk(bedfile, bimfile, famfile, *pM);
  return pM;
}
