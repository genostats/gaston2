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
Rcpp::IntegerVector writeBedFileDisk_input(std::string path, size_t n_ind, size_t n_snp, size_t index) {
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


// [[Rcpp::export]]
Rcpp::IntegerVector writeBedFileDisk_newfile(std::string path, size_t n_ind, size_t n_snp) {
  // TODO: check if file exists, refuse to erase
  std::ofstream file_test(path, std::ifstream::out); //will erase its content cos not opened in append
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  // Adding magic numbers for mio to see that file is not empty
  file_test.put(108);
  file_test.put(27);
  file_test.put(1);
  file_test.flush();
  std::error_code error;
  mio::mmap_sink file_ = mio::make_mmap_sink(path, 0, mio::map_entire_file, error);
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
  }

  // Fails if file is empty, so magic byte put first
  std::shared_ptr<mio::mmap_sink> file_ptr = std::make_shared<mio::mmap_sink>(std::move(file_));
  // char *data = file_ptr->data();
  // // check magic number
  // char magic[3] = {108, 27, 1};
  // for (int i = 0; i < 3; i++)
  // {
  //   // std::cout << int(*data) << ' '; // this prints 108 27 1
  //   *data++ = magic[i] ;
  // }

  // creating the vector of SNP to write, will be zero
  size_t bytes_for_SNP = n_ind/4 + ((n_ind%4 == 0u)?0:1);
  std::vector<uint8_t> SNP2write(bytes_for_SNP);
  for (int i = 0; i < bytes_for_SNP; i++) {
    SNP2write[i] = 50; 
  }

  // For clarity :
  std::cout << "This is the SNP to write in file ";
  for (uint8_t i: SNP2write)
    std::cout << i;
  std::cout <<'\n';

  for (size_t i = 0; i < n_snp; i++) {
    try {
      SNPWriteDisk(n_ind, file_ptr, i, SNP2write, path);
    } catch(const std::exception& e) {
      std::cerr << e.what() << '\n';
    }
  }
  
  // now reading
  SNPmatrix M = readBedFileDisk(path, n_ind, n_snp);
  std::vector<int> res;
  for(auto v : M.SNPs) {
    res.push_back(v->sum());
  }
  return Rcpp::wrap(res);
}


// [[Rcpp::export]]
void testWNF(std::string path) {
  // TODO: check if file exists, refuse to erase
  std::ofstream file_test(path, std::ifstream::out); //will erase its content cos not opened in append
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  // Adding magic numbers for mio to see that file is not empty
  file_test.put(108);
  file_test.put(27);
  file_test.put(1);
  file_test.flush();
  std::error_code error;
  mio::mmap_sink file_ = mio::make_mmap_sink(path, 0, mio::map_entire_file, error);
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
  }

  char * DD = file_.data();

  DD[0] = 1;
  resizing_stream(file_test, 5);
  DD[3] = 8;
  DD[4] = 16;
  
  resizing_stream(file_test, 5);
  DD[8] = 32;
  DD[9] = 24;

  file_test.close();
}

// Rcpp::IntegerVector writeBedFileDisk_fail(std::string path, size_t n_ind, size_t n_snp, size_t index) {
//   // first checking if works with a file not opened


// }