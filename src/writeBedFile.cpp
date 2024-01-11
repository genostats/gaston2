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
#include "SNPwriteDisk.h"


/*
SNPmatrix readBedFileDisk(std::string path, size_t n_ind, size_t n_snp, size_t index) {
  std::ifstream file_test(path, std::ifstream::binary);
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  std::error_code error;
  mio::mmap_sink * file_ = mio::make_mmap_sink(path, 0, mio::map_entire_file, error);
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
  }

  // check magic number
  char magic[3];
  for (int i = 0; i < 3; i++)
  {
    std::cout << int(*data) << ' ';
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
} */