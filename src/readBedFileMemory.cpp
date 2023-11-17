#include "SNPmatrix.h"
#include "SNPvectorMemory.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <memory>

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
    std::shared_ptr<SNPvectorMemory> snpVec(new SNPvectorMemory(n_ind));
    size_t n = snpVec->nbChars();
    uint8_t * data = snpVec->data();
    file.read(reinterpret_cast<char *>(data), n);
    M.push_back(snpVec);
  }
  
  file.close();
  return M;
}

