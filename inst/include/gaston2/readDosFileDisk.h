/**
 * @file readDosFileDisk.cpp
 * @author Hervé Perdry
 * @brief Main file with functions to interact with SNPvector and SNPmatrix
 * @date 2023-12-11
 * 
 */
 #ifndef READOSFILEDISK_H
 #define READOSFILEDISK_H

#include "SNPmatrix.h"
#include "SNPdosageDisk.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <memory>
#include "mio.hpp"
#include <cstring>
#include <Rcpp.h>

/**
 * @brief Reading a .dosf file, storing SNPs in a SNPmatrix returned
 * 
 * @param dosfile The name of the dos file
 * @param bimfile The name of the bim file
 * @param bamfile The name of the fam file
 * @param M reference to an originally empty SNPmatrix
 */

inline void readDosFileDisk(std::string dosfile, std::string bimfile, std::string famfile, SNPmatrix<SNPdosage> & M) {
  // first read bim and fam file, to determine size
  //std::cout << "bim file is " << bimfile;
  
  M.readFamFile(famfile);
  M.readBimFile(bimfile);

  size_t n_snp = M.getSNPStats().nrow();
  size_t n_ind = M.getIndStats().nrow();

  // now ready to read dos file
  std::ifstream file_test(dosfile, std::ifstream::binary);
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  std::error_code error;
  mio::mmap_source file_ = mio::make_mmap_source(dosfile, 0, mio::map_entire_file, error);
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
  }
  std::shared_ptr<mio::mmap_source> file_ptr = std::make_shared<mio::mmap_source>(std::move(file_)); // don't know if good idea, creates a NEW sink
  // DON'T USE FILE_ FROM NOW ON COS IT'S NULL, WAS MOVED
  const float* data = reinterpret_cast<const float*>(file_ptr->data());// const necessary bcos read only 
  

  for(size_t i = 0; i < n_snp; i++) {
    //makes a shared_ptr on a vector of snips 
    std::shared_ptr<SNPdosageDisk<mio::access_mode::read>> snpVec(new SNPdosageDisk<mio::access_mode::read>(n_ind, file_ptr, i));
    M.push_back(snpVec);
  }
}

inline SNPmatrix<SNPdosage> readDosFileDisk(std::string bedfile, std::string bimfile, std::string famfile) {
  SNPmatrix<SNPdosage> M;
  readDosFileDisk(bedfile, bimfile, famfile, M);
  return M;
}

#endif /*READOSFILE_H*/