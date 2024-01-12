#include "SNPmatrix.h"
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
#include <Rcpp.h>
#include "readBedFileMemory.h"
#include "SNPwriteDisk.h"


// [[Rcpp::export]]
Rcpp::IntegerVector writeBedFileDisk(std::string path, size_t n_ind, size_t n_snp, size_t index) {
  std::ifstream file_test(path, std::ifstream::binary);
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  std::error_code error;
  mio::mmap_sink file_ = mio::make_mmap_sink(path, 0, mio::map_entire_file, error);
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
  }

  // HERE ADD shared pointer to pass on 
  // should I deal with it with a wait ? if SNP by SNP should be good
  std::shared_ptr<mio::mmap_sink> file_ptr = std::make_shared<mio::mmap_sink>(std::move(file_));
  char *data = file_ptr->data();
  // check magic number
  char magic[3];
  for (int i = 0; i < 3; i++)
  {
    // std::cout << int(*data) << ' '; // this prints 108 27 1
    magic[i] = *data++;
  }
  if(magic[0] != 108 || magic[1] != 27) {
    throw std::runtime_error("Not a bed file");
  }
  if(magic[2] != 1) {
    throw std::runtime_error("Not a bed file in SNP major mode");
  }

  std::vector<uint8_t> SNP2write;
  Rcpp::Function readline("readline");

  size_t bytes_for_SNP = n_ind/4 + ((n_ind%4 == 0u)?0:1);
  
  Rcpp::Rcout << "Enter a SNP of size (so ... byte)" << bytes_for_SNP;
  std::string input = Rcpp::as<std::string>(readline("> "));

  uint8_t value;
  int l = input.length();
  if ( l != bytes_for_SNP) throw std::runtime_error("The input isn't of the right size !\n");
  for (int i = 0; i < l; i++) {
    value = input.at(i);
    SNP2write.push_back(value - 48); // 48 = '0'
  }

  // For clarity :
  std::cout << "This is the SNP to write in file ";
  for (char i: SNP2write)
    std::cout << i;
  std::cout <<'\n';

  SNPWriteDisk(n_ind, file_ptr, index, SNP2write, path);

  // now reading
  SNPmatrix M = readBedFileDisk(path, n_ind, n_snp);
  std::vector<int> res;
  for(auto v : M.SNPs) {
    res.push_back(v->sum());
  }
  return Rcpp::wrap(res);
}