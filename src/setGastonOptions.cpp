#include "gastonOptions.h"
#include <Rcpp.h>
#include "debug.h"

// [[Rcpp::export]] 
void setGastonOptions(Rcpp::List L) {

  Rcpp::IntegerVector autosomes_ = L["autosomes"];
  std::set<int> autosomes = std::set<int>(autosomes_.begin(), autosomes_.end());

  Rcpp::IntegerVector x_ = L["x"];
  std::set<int> x = std::set<int>(x_.begin(), x_.end());

  Rcpp::IntegerVector y_ = L["y"];
  std::set<int> y = std::set<int>(y_.begin(), y_.end());

  Rcpp::IntegerVector mt_ = L["mt"];
  std::set<int> mt = std::set<int>(mt_.begin(), mt_.end());

  setGastonOptions(autosomes, x, y, mt);
}
