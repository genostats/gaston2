#include <Rcpp.h>
#include "SNPmatrix.h"
#include "SNPvector.h"
#include "readBedFileMemory.h"
#include "debug.h"
#include <vector>

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
