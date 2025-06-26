#include "SNPmatrix.h"
#include "SNPvectorDisk.h"
#include <iostream>
#include <memory> // shared_ptr
#include <cstring>
#include <Rcpp.h>

#ifndef _BINDINDMATRIX_DISK_
#define _BINDINDMATRIX_DISK_


// filling up a new SNPmatrix on disk (in the file specified by path_str) with every individuals from the first and second matrix
void bindIndstoSNPmatrixDisk(const SNPmatrix<SNPvector> &first, const SNPmatrix<SNPvector> &second, std::string path_str, SNPmatrix<SNPvectorDisk<mio::access_mode::write>> &newMat) {

  std::error_code error;

  if (first.size() != second.size())
  throw std::logic_error("You cannot merge 2 SNPmatrix with mismatched number of SNPs (maybe later could add NAs)");


  /* FIRST : check if file exists, if it does, abort to not overwrite */
  const char * path_c = path_str.c_str();
  FILE *check = fopen(path_c, "rb");
  if (check) {
      fclose(check);
      throw std::runtime_error("The new matrix will overwrite an existing file, Aborting !");
  }

  /* THEN if it doesn't exist, go prepare it for mio with .bed's magic bytes */
  FILE *f = fopen(path_c, "wb");
  if (!f)
  {
    throw std::runtime_error("Failed to open file for writing");
  }
  // Adding magic numbers to Identify a bed file in SNP major mode
  fputc(108, f);
  fputc(27, f);
  fputc(1, f);


  const std::vector<std::shared_ptr<SNPvector>> firstSNPs = first.getSNPs();
  const std::vector<std::shared_ptr<SNPvector>> secondSNPs = second.getSNPs();

  const int total_nb_inds = first.nbInds() + second.nbInds();
  int total_nb_char = total_nb_inds / 4 + ((total_nb_inds % 4 == 0)? 0 : 1);

  // (+ 3 for the 3 magic bytes)
  int to_add = total_nb_char * first.size() + 3 ;
  if (fseek(f, to_add - 1, SEEK_SET) != 0)
  {
    fclose(f);
    throw std::runtime_error("Error when resizing file");
  }
  fputc(0, f); // this is what will size it up
  fclose(f);

  /* NOW : mapping the created file with mio, */
  std::shared_ptr<mio::mmap_sink> file_ptr = std::make_shared<mio::mmap_sink>(mio::make_mmap_sink(path_str, 0, mio::map_entire_file, error));
  if (error)
  {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg);
  }

  // maybe could also do a constructor or a function on the SNPmatrix level ?
  for (size_t i = 0; i < firstSNPs.size(); i++) {
    newMat.push_back(std::make_shared<SNPvectorDisk<mio::access_mode::write>>(firstSNPs[i], secondSNPs[i], file_ptr, i));
  }

  // Individuals stats are still good (N0, N1, N2, NAs also)
  // so I can just fuse them
  newMat.setIndStats(DataStruct(first.getIndStats(), second.getIndStats()));
  // SNP stats need to be recomputed
  newMat.setSnpStats(first.getSNPStats());
  newMat.setMode(first.getMode());

}

SNPmatrix<SNPvectorDisk<mio::access_mode::write>> bindIndstoSNPmatrixDisk(const SNPmatrix<SNPvector> &first, const SNPmatrix<SNPvector> &second, std::string path) {
  SNPmatrix<SNPvectorDisk<mio::access_mode::write>> M;
  bindIndstoSNPmatrixDisk(first, second, path, M);
  return M;
}

#endif 
