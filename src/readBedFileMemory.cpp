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


/**
 * @brief Reading a bed file, storing SNPs in a SNPmatrix returned
 * 
 * @param filename The name of the file to be opened
 * @param n_ind The number of individuals/samples
 * @param n_snp The number of SNP to read from the file and to load into the Matrix.
 * 
 * @return SNPmatrix, stocking shared_ptrs to SNPvectorMemory 
 */
SNPmatrix readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp) {
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
    std::shared_ptr<SNPvectorMemory> snpVec(new SNPvectorMemory(n_ind));
    size_t n = snpVec->nbChars();
    uint8_t * data = snpVec->data();
    file.read(reinterpret_cast<char *>(data), n);
    M.push_back(snpVec);
  }
  
  file.close();
  return M;
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
SNPmatrix readBedFileDisk(std::string path, size_t n_ind, size_t n_snp) {
  std::ifstream file_test(path, std::ifstream::binary);
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  std::error_code error;
  mio::mmap_source file_ = mio::make_mmap_source(path, 0, mio::map_entire_file, error);
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
    // std::cout << int(*data) << ' ';
    magic[i] = *data++;
  }
  std::cout << '\n';
  if(magic[0] != 108 || magic[1] != 27) {
    throw std::runtime_error("Not a bed file");
  }
  if(magic[2] != 1) {
    throw std::runtime_error("Not a bed file in SNP major mode");
  }
  
  SNPmatrix M;
  auto file_offset = 3; // BCOS MAGIC BYTES 
  for(size_t i = 0; i < n_snp; i++) {
    std::shared_ptr<SNPVectorDisk> snpVec(new SNPVectorDisk(n_ind,file_ptr, i));
    size_t n = snpVec->nbChars(); // func inherited from SNPVec, gives back sizof vec

    // TODO : check if this part works
    // maybe no need to give it, already done in constructor...


    // uint8_t * data = snpVec->data(); // this data is ptr to first char of SNPVector
    // //Copies count n bytes from file_ptr->data to data, vector of SNPVec. Both are reinterpreted as arrays of unsigned char. 
    // std::memcpy(data, file_ptr->data() + file_offset, n);

    M.push_back(snpVec);
    // // Increment the file offset based on the size of the data
    file_offset += n;
  }
  return M;
}