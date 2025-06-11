#include <Rcpp.h>

#include <iostream>

#include "SNPmatrix.h"
#include "SNPvectorDisk.h"  //should include also SNPvector

/* Function that eats a SNPmatrix<SNPvector>
and fills back an instance of SNPmatrix<SNPvectorDisk> specifically
with all the data and metadata from the original matrix.
Also creating a corresponding file. /!\ This matrix HAS to be with SNPvectorDisk<mio::access_mode::write>
*/
void ToSNPmatrixDisk(const SNPmatrix<SNPvector> &other, std::string newfile_name, SNPmatrix<SNPvectorDisk<mio::access_mode::write>> &newMat) {
  std::error_code error;

  /* FIRST : check if file exists, if it does, abort to not overwrite */
  const char *path_c = newfile_name.c_str();
  FILE *check = fopen(path_c, "rb");
  if (check) {
    fclose(check);
    throw std::runtime_error("The new matrix will overwrite an existing file, Aborting !");
  }

  /* THEN if it doesn't exist, go prepare it for mio with .bed's magic bytes */
  FILE *f = fopen(path_c, "wb");
  if (!f) {
    throw std::runtime_error("Failed to open file for writing");
  }
  // Adding magic numbers to Identify a bed file in SNP major mode
  fputc(108, f);
  fputc(27, f);
  fputc(1, f);

  const std::vector<std::shared_ptr<SNPvector>> otherSNPs = other.getSNPs();

  // (+ 3 for the 3 magic bytes)
  int to_add = (other.nbInds() / 4 + ((other.nbInds() % 4 == 0) ? 0 : 1)) * otherSNPs.size() + 3;
  if (fseek(f, to_add - 1, SEEK_SET) != 0) {
    fclose(f);
    throw std::runtime_error("Error when resizing file");
  }
  fputc(0, f);  // this is what will size it up
  fclose(f);

  /* NOW : mapping the created file with mio, */
  std::shared_ptr<mio::mmap_sink> file_ptr = std::make_shared<mio::mmap_sink>(mio::make_mmap_sink(newfile_name, 0, mio::map_entire_file, error));
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg);
  }

  size_t SNPindex = 0;

  for (const auto &snp : otherSNPs) {
    newMat.push_back(std::make_shared<SNPvectorDisk<mio::access_mode::write>>(snp, file_ptr, SNPindex++));
  }
  // extract stats now and set stats_set_ to true
  newMat.setIndStats(other.getIndStats());
  newMat.setSnpStats(other.getSNPStats());
}

SNPmatrix<SNPvectorDisk<mio::access_mode::write>> ToSNPmatrixDisk(SNPmatrix<> other, std::string newfile_name) {
    SNPmatrix<SNPvectorDisk<mio::access_mode::write>> M;
    ToSNPmatrixDisk(other, newfile_name, M);
    return M;
  }

// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPvectorDisk<mio::access_mode::write>>> ToSNPmatrixDisk_(Rcpp::XPtr<SNPmatrix<SNPvector>> pM, std::string newfile_name) {
  Rcpp::XPtr<SNPmatrix<SNPvectorDisk<mio::access_mode::write>>> pnewMat(new SNPmatrix<SNPvectorDisk<mio::access_mode::write>>);
  ToSNPmatrixDisk(*pM, newfile_name, *pnewMat);
  return pnewMat;
}