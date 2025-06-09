#include "SNPmatrix.h"
#include <iostream>
#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::IntegerVector dimDosematrix(Rcpp::XPtr<SNPmatrix<>> pM) {
  return Rcpp::IntegerVector::create(pM->nbInds(), pM->nbSNPs());
}