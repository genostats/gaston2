#include "SNPmatrix.h"
#include "SNPvectorMemory.h"
#include "bindIndstoSNPmatrixMemory.h"
#include <Rcpp.h>

//R exported function to have a new SNPmatrix<SNPvector> of combined inds 
// from first and second matrix
// could be recalled rbind !!
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPvectorMemory>> bindIndstoSNPmatrixMemory_(Rcpp::XPtr<SNPmatrix<SNPvector>> first, Rcpp::XPtr<SNPmatrix<SNPvector>> second) {
  Rcpp::XPtr<SNPmatrix<SNPvectorMemory>> pnewMat(new SNPmatrix<SNPvectorMemory>);
  bindIndstoSNPmatrixMemory(*first, *second, *pnewMat);
  return pnewMat;
}