#include "SNPmatrix.h"
#include "hwe.h"
#include "set_hwe.h"
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
void set_hwe(Rcpp::XPtr<SNPmatrix<>> pM, int test = 0, double usefloat = false) {

  // TODO settle this. SNPStats always in double?
  // we really need exportSNPStats and not only computeSNPStats because of chrX...
  pM->exportSNPStats<double>(false); 

  DataStruct & ds = pM->getSNPStats();

  if(usefloat) {
    switch(test) {
      case 0:
        set_hwe<HWEtest::chisquare, float>(ds); 
        break;
      case 1:
        set_hwe<HWEtest::chisquare_yates, float>(ds); 
        break;
      case 2:
        set_hwe<HWEtest::exact, float>(ds); 
        break;
      default:
        Rcpp::stop("unknown test code");
     }
  } else {
    switch(test) {
      case 0:
        set_hwe<HWEtest::chisquare, double>(ds); 
        break;
      case 1:
        set_hwe<HWEtest::chisquare_yates, double>(ds); 
        break;
      case 2:
        set_hwe<HWEtest::exact, double>(ds); 
        break;
      default:
        Rcpp::stop("unknown test code");
    }
  }
}
