#include "isAutosome.h"
#include <Rcpp.h>
#include "debug.h"

// [[Rcpp::export]]
Rcpp::LogicalVector isAutosome_(Rcpp::IntegerVector chr) {
  return Rcpp::wrap( isAutosome(chr) );
}

// [[Rcpp::export]]
Rcpp::LogicalVector isX_(Rcpp::IntegerVector chr) {
  return Rcpp::wrap( isX(chr) );
}

// [[Rcpp::export]]
Rcpp::LogicalVector isY_(Rcpp::IntegerVector chr) {
  return Rcpp::wrap( isY(chr) );
}

// [[Rcpp::export]]
Rcpp::LogicalVector isMt_(Rcpp::IntegerVector chr) {
  return Rcpp::wrap( isMt(chr) );
}

