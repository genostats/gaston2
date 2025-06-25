#include "SNPmatrix.h"
#include "SNPdosageMemory.h"
#include "bindIndstoDosagematrixMemory.h"
#include <Rcpp.h>

//R exported function to have a new SNPmatrix<SNPdosageMemory> of combined inds 
// from first and second matrix
// could be recalled rbind !!
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosageMemory>> bindIndstoDosagematrixMemory_(Rcpp::XPtr<SNPmatrix<SNPdosage>> first, Rcpp::XPtr<SNPmatrix<SNPdosage>> second) {
  Rcpp::XPtr<SNPmatrix<SNPdosageMemory>> pnewMat(new SNPmatrix<SNPdosageMemory>);
  bindIndstoDosagematrixMemory(*first, *second, *pnewMat);
  return pnewMat;
}