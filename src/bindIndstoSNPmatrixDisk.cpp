#include "SNPmatrix.h"
#include "SNPvectorDisk.h"
#include "bindIndstoSNPmatrixDisk.h"
#include <Rcpp.h>

//R exported function to have a new SNPmatrix<SNPvector> of combined inds 
// from first and second matrix
// could be recalled rbind !!
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPvectorDisk<mio::access_mode::write>>> bindIndstoSNPmatrixDisk_(Rcpp::XPtr<SNPmatrix<SNPvector>> first, Rcpp::XPtr<SNPmatrix<SNPvector>> second, std::string newmatrix_path) {
  Rcpp::XPtr<SNPmatrix<SNPvectorDisk<mio::access_mode::write>>> pnewMat(new SNPmatrix<SNPvectorDisk<mio::access_mode::write>>);
  bindIndstoSNPmatrixDisk(*first, *second, newmatrix_path, *pnewMat);
  return pnewMat;
}