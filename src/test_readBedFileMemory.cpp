#include <Rcpp.h>
#include "SNPmatrix.h"
#include "SNPvector.h"
#include "readBedFileMemory.h"
#include "debug.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector test_readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp) {
  SNPmatrix M = readBedFileMemory(filename, n_ind, n_snp);
  std::vector<int> res;
  for(auto v : M.SNPs) {
    res.push_back(v->sum());
  }
  return wrap(res);
}

// [[Rcpp::export]]
IntegerVector test_readBedFileDisk(std::string filename, size_t n_ind, size_t n_snp) {
  SNPmatrix M = readBedFileDisk(filename, n_ind, n_snp);
  std::vector<int> res;
  for(auto v : M.SNPs) {
    res.push_back(v->sum());
  }
  return wrap(res);
}

// [[Rcpp::export]]
IntegerVector test_delete(std::string filename, size_t n_ind, size_t n_snp) {
  SNPmatrix M = readBedFileDisk(filename, n_ind, n_snp);
  std::vector<int> res;
  for(auto v : M.SNPs) {
    res.push_back(v->sum());
  }
  SNPmatrix MM = readBedFileMemory(filename, n_ind, n_snp);
  for(auto v : MM.SNPs) {
    M.push_back(v);
  }
  std::cout << "Is the 3rd SNP on disk ? (yes) : " << M.onDisk(3);
  std::cout << "\nAnd this one ? (no) : " << M.onDisk(n_snp + 3) << "\n";
  while (M.size() != 0) {
    M.deleteSNP();
  }
  std::cout << "Is the file destroyed before or after ?\n";
  return wrap(res);
}

// test reading only first bit of SNPs in the same file, opened in 2 different ways
// [[Rcpp::export]]
IntegerVector test_read_parts(std::string filename, size_t n_ind, size_t n_snp) {
  SNPmatrix Mem = readBedFileMemory(filename, n_ind, n_snp);
  std::vector<int> res;
  for (auto v : Mem.SNPs) {
    res.push_back(v->data()[0]);
  }
  for (int i = 0; i < 10; i++) {
    res.push_back(0);
  }
  SNPmatrix Dis = readBedFileDisk(filename, n_ind, n_snp);
  for (auto v : Dis.SNPs) {
    res.push_back(*v->data());
  }
  return wrap(res);
}