#include "gastonOptions.h"
#include <Rcpp.h>
#include "debug.h"

// [[Rcpp::export]]
Rcpp::List getGastonOptions_() {
  gastonOptions & opt = getGastonOptions();

  Rcpp::List L;
  L["autosomes"] = Rcpp::wrap(std::vector<int>(opt.autosomes.begin(), opt.autosomes.end())); 
  L["x"] = Rcpp::wrap(std::vector<int>(opt.x.begin(), opt.x.end())); 
  L["y"] = Rcpp::wrap(std::vector<int>(opt.y.begin(), opt.y.end())); 
  L["mt"] = Rcpp::wrap(std::vector<int>(opt.mt.begin(), opt.mt.end())); 
  
  return L;
}

