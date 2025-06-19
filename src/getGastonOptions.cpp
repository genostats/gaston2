#include "gastonOptions.h"
#include <Rcpp.h>
#include "debug.h"

// [[Rcpp::export]]
Rcpp::List getGastonOptions_() {
  gastonOptions & opt = getGastonOptions();

  Rcpp::List L;
  L["autosomes"] = Rcpp::wrap(opt.autosomes); 
  L["x"] = Rcpp::wrap(opt.x); 
  L["y"] = Rcpp::wrap(opt.y); 
  L["mt"] = Rcpp::wrap(opt.mt); 
  
  return L;
}

