#include "SNPmatrix.h"
#include "SNPdosageDisk.h"
#include "bindIndstoDosagematrixDisk.h"
#include <Rcpp.h>

//R exported function to have a new SNPmatrix<SNPdosageDisk> of combined inds 
// from first and second matrix
// could be recalled rbind !!
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosageDisk<mio::access_mode::write>>> bindIndstoDosagematrixDisk_(Rcpp::XPtr<SNPmatrix<SNPdosage>> first, Rcpp::XPtr<SNPmatrix<SNPdosage>> second, std::string path) {
  Rcpp::XPtr<SNPmatrix<SNPdosageDisk<mio::access_mode::write>>> pnewMat(new SNPmatrix<SNPdosageDisk<mio::access_mode::write>>);
  bindIndstoDosagematrixDisk(*first, *second, path, *pnewMat);
  return pnewMat;
}