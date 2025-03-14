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
  std::vector<unsigned int> res;
  for(auto v : M.SNPs) {
    res.push_back(v->sum());
  }
  return wrap(res);
}

// [[Rcpp::export]]
IntegerVector test_readBedFileDisk(std::string filename, size_t n_ind, size_t n_snp) {
  SNPmatrix M = readBedFileDisk(filename, n_ind, n_snp);
  std::vector<unsigned int> res;
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

IntegerVector loop_sum(SNPmatrix matrix) {
  std::vector<int> res;
  for(auto v : matrix.SNPs) {
    //== sum(503)
    res.push_back(v->sum());
  }
  return wrap(res);
}


// [[Rcpp::export]]
IntegerVector test_readModes(std::string filename, size_t n_ind, size_t n_snp) {
  std::cout << " reading : " << filename << "\n n_ind : " << n_ind << "\n n_snp : " << n_snp << "\n";
  std::cout << " reading in Numeric Mode \n";
  uint8_t Num[4] = {0, 1, 2, 3};
  SNPmatrix M = readBedFileMemory(filename, n_ind, n_snp, Num);
  IntegerVector res = loop_sum(M);
  std::cout << " reading in Centered (not implemented yet) Mode \n";
  uint8_t Cent[4] = { 0, 0, 0, 0};
  M = readBedFileMemory(filename, n_ind, n_snp, Cent);
  IntegerVector res2 = loop_sum(M);
  for (auto i : res2) res.push_back(i);
  std::cout << " reading in Standardized (not implemented yet) Mode \n";
  uint8_t Std[4] = {0, 0, 0, 0};
  M = readBedFileMemory(filename, n_ind, n_snp, Std);
  IntegerVector res3 = loop_sum(M);
  for (auto i : res3) res.push_back(i);
  std::cout << " reading in PLINK Mode \n";
  M = readBedFileMemory(filename, n_ind, n_snp);
  IntegerVector res4 = loop_sum(M);
  for (auto i : res4) res.push_back(i);
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
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503) { // cas où je dépasse SNP[0]
    nbSNPs = (n/503) + 1;
    n = n%503;
  } //parce que 503 individus dans le file hardcodé
  //std::cout << "Nb of SNPs: " << nbSNPs << " and nb of bits to read on last :" << n <<"\n";
  SNPmatrix M = readBedFileMemory(file_hardcode, (nbSNPs > 1)? 503 : n, nbSNPs);
  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i]; //peut pas déréférencer là parce qu'instancie la classe abstraite SNPvector
    if (i == nbSNPs -1) {
      S += vec->sum(n);
      //std::cout << "final S: " << S <<"\n";
      return S;
    }
    S += vec->sum();
    //std::cout << "S for now: " << S <<"\n";

  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_1(unsigned int n) {
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503) { // cas où je dépasse SNP[0]
    nbSNPs = (n/503) + 1;
    n = n%503;
  }
  SNPmatrix M = readBedFileMemory(file_hardcode, (nbSNPs > 1)? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i];
    if (i == nbSNPs -1) {
      S += vec->sum(n);
      return S;
    }
    for(unsigned int a : *vec)
      S += a;
  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_2(unsigned int n) {
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503) {
    nbSNPs = (n/503) + 1;
    n = n%503;
  }
  SNPmatrix M = readBedFileMemory(file_hardcode, (nbSNPs > 1)? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i];
    if (i == nbSNPs -1) {
      S += vec->sum(n);
      return S;
    }
  for(auto pa = vec->begin(); pa != vec->end(); ++pa)
      S += *pa;
  }
  return -1;
}

// ON DISK NOW, mean = 15 microseconds

// [[Rcpp::export]]
unsigned int test_performance_iterator_disk(unsigned int n) {
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503) { // cas où je dépasse SNP[0]
    nbSNPs = (n/503) + 1;
    n = n%503;
  }
  SNPmatrix M = readBedFileDisk(file_hardcode, (nbSNPs > 1)? 503 : n, nbSNPs);
  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i];
    if (i == nbSNPs -1) {
      S += vec->sum(n);
      return S;
    }
    S += vec->sum();
  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_1d(unsigned int n) {
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503) { // cas où je dépasse SNP[0]
    nbSNPs = (n/503) + 1;
    n = n%503;
  }
  SNPmatrix M = readBedFileMemory(file_hardcode, (nbSNPs > 1)? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i];
    if (i == nbSNPs -1) {
      S += vec->sum(n);
      return S;
    }
    for(unsigned int a : *vec)
      S += a;
  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_2d(unsigned int n) {
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503) {
    nbSNPs = (n/503) + 1;
    n = n%503;
  }
  SNPmatrix M = readBedFileDisk(file_hardcode, (nbSNPs > 1)? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i];
    if (i == nbSNPs -1) {
      S += vec->sum(n);
      return S;
    }
  for(auto pa = vec->begin(); pa != vec->end(); ++pa)
      S += *pa;
  }
  return -1;
}

// TODO : 

// [[Rcpp::export]]
IntegerVector test_snp_stats(unsigned int n) {
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503) {
    nbSNPs = (n/503) + 1;
    n = n%503;
  }
  SNPmatrix M = readBedFileMemory(file_hardcode, (nbSNPs > 1)? 503 : n, nbSNPs);

  std::vector<unsigned int> res(4, 0);  // Initialize sum vector {0, 0, 0, 0}

  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i]; 
    auto stats = vec->compute_stats();

    for (int j = 0; j < 4; j++) {
      res[j] += stats[j];  // Element-wise sum
    }
    // FOR now no way to stop on a part of the SNP
    if (i == nbSNPs -1) {
      return wrap(res);
    }

  }
  return -1;
}

// [[Rcpp::export]]
IntegerVector test_snp_stats_d(unsigned int n) {
  if (n > 503) n = 503;
  SNPmatrix M = readBedFileDisk(file_hardcode, n, 1);
  auto vec = M.SNPs[0];
  auto stats = vec->compute_stats();
  std::vector<unsigned int> res(stats, stats + 4);
  return wrap(res);
}

/*******************************
 *    Test SNP sum comparison  *
 ********************************/

// [[Rcpp::export]]
IntegerVector test_sums(unsigned int n){

  unsigned int sum_default = test_performance_iterator_default(n);
  auto sum_1 = test_performance_iterator_1(n);
  auto sum_2 = test_performance_iterator_2(n);
  auto sum_disk = test_performance_iterator_disk(n);
  auto sum_1d = test_performance_iterator_1d(n);
  auto sum_2d = test_performance_iterator_2d(n);

  if (std::set<unsigned int>{sum_default, sum_1, sum_2, sum_disk, sum_1d, sum_2d}.size() == 1) { // set rm duplicates, if size == 1 only 1 value
    std::cout << "All sums are equal!\n";
  }

  IntegerVector stats = test_snp_stats(n);
  return stats;
}


/*******************************
 *    Test mu and sigma  *
 ********************************/

// [[Rcpp::export]]
NumericVector test_mu_sigma(unsigned int n) {
  int nbSNPs = n;

  SNPmatrix M = readBedFileMemory(file_hardcode, 503, nbSNPs);

  std::vector<double> res;

  for (int i = 0; i < nbSNPs; i++) {
    auto vec = M.SNPs[i]; 
    // NEED TO COMPUTE BEFORE 
    auto stats = vec->compute_stats();
    //std::cout << "This is stats for SNP[" << i << "] : " << stats[0] << ", " << stats[1] << ", " << stats[2];

    auto musi = vec->compute_mu_sigma();
    double sigma = vec->sigma();
    double mu = vec->mu();
    res.push_back(mu);
    res.push_back(sigma);
    res.push_back(0);
    std::cout << "These are stats for SNP[" << i << "], mu: " << mu << ", sigma: " << sigma << "\n";
  }
  return wrap(res);
}
