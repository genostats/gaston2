#include "gastonOptions.h"
#include <Rcpp.h>
#include "debug.h"

// [[Rcpp::export]] 
void setGastonOptions(Rcpp::List L) {
  std::vector<int> autosomes = Rcpp::as<std::vector<int>>(L["autosomes"]);
  std::vector<int> x = Rcpp::as<std::vector<int>>(L["x"]);
  std::vector<int> y = Rcpp::as<std::vector<int>>(L["y"]);
  std::vector<int> mt = Rcpp::as<std::vector<int>>(L["mt"]);
  setGastonOptions(autosomes, x, y, mt);
}
