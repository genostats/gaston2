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

//when testing using : "extdata/LCT.bed" n_ind : 503 n_snp : 607

// [[Rcpp::export]]
IntegerVector test_readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp) {
  std::cout << " reading : " << filename << "\n n_ind : " << n_ind << "\n n_snp : " << n_snp << "\n";
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

/************************
 *    Test SNP reading  *
 ************************/


// TODO : to test with numerous SNPs, so other rows in the bed matrix

const char file_hardcode[68] = "/home/ju/R/x86_64-pc-linux-gnu-library/4.4/snipsnop/extdata/LCT.bed";

// ON MEMORY, mean = 4.5 microseconds on average

// [[Rcpp::export]]
unsigned int test_performance_iterator_default(unsigned int n) {
  if (n > 503) n = 503; //parce que 503 individus dans le file hardcodé
  SNPmatrix M = readBedFileMemory(file_hardcode, n, 1);
  auto vec = M.SNPs[0]; //peut pas déréférencer là parce qu'instancie la classe abstraite SNPvector
  unsigned int S = vec->sum();
  return S;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_1(unsigned int n) {
  if (n > 503) n = 503;
  SNPmatrix M = readBedFileMemory(file_hardcode, n, 1);
  auto vec = M.SNPs[0];
  unsigned int S = 0;
  for(unsigned int a : *vec)
    S += a;
  return S;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_2(unsigned int n) {
  if (n > 503) n = 503;
  SNPmatrix M = readBedFileMemory(file_hardcode, n, 1);
  auto vec = M.SNPs[0];
  unsigned int S = 0;
  // TODO : check with type SNPvectorMemory::Iterator pa ???
  for(auto pa = vec->begin(); pa != vec->end(); ++pa)
    S += *pa;
  return S;
}

// ON DISK NOW, mean = 15 microseconds

// [[Rcpp::export]]
unsigned int test_performance_iterator_disk(unsigned int n) {
  if (n > 503) n = 503; //parce que 503 individus dans le file hardcodé
  SNPmatrix M = readBedFileDisk(file_hardcode, n, 1);
  auto vec = M.SNPs[0]; //peut pas déréférencer là parce qu'instancie la classe abstraite SNPvector
  unsigned int S = vec->sum();
  return S;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_1d(unsigned int n) {
  if (n > 503) n = 503;
  SNPmatrix M = readBedFileDisk(file_hardcode, n, 1);
  auto vec = M.SNPs[0];
  unsigned int S = 0;
  for(unsigned int a : *vec)
    S += a;
  return S;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_2d(unsigned int n) {
  if (n > 503) n = 503;
  SNPmatrix M = readBedFileDisk(file_hardcode, n, 1);
  auto vec = M.SNPs[0];
  unsigned int S = 0;
  // TODO : check with type SNPvectorMemory::Iterator pa ???
  for(auto pa = vec->begin(); pa != vec->end(); ++pa)
    S += *pa;
  return S;
}