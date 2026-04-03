#include "SNPmatrix.h"
#include "SNPvectorDisk.h"
#include <iostream>
#include <memory> // shared_ptr
#include "mio.hpp"
#include <cstring>

#ifndef _EXTRACTINDMATRIX_DISK_
#define _EXTRACTINDMATRIX_DISK_

template <typename SNPvectorClass, typename intVec>
void extractIndsfromSNPmatrixDisk(const SNPmatrix<SNPvectorClass> &other, const intVec &keep, std::string path_str, 
                                  SNPmatrix<SNPvectorClass> & newMat) {

  std::error_code error;

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

  const std::vector<std::shared_ptr<SNPvector>> otherSNPs = other.getSNPs();

  // (+ 3 for the 3 magic bytes)
  int to_add = (keep.size() / 4 + ((keep.size() % 4 == 0u) ? 0 : 1)) * otherSNPs.size() + 3 ;
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

  size_t SNPindex = 0;

  for (const auto &snp : otherSNPs){
    // not populating stats now
    newMat.push_back(std::make_shared<SNPvectorDisk<mio::access_mode::write>>(snp, file_ptr, SNPindex++, keep));
  }
  //extract stats now and set indStatsComputed to true
  DataStruct original_dt = other.getIndStats();
  newMat.setIndStats(DataStruct(original_dt, keep));
  // and set them to be as "complete" as the mother matrix
  newMat.setindStatscomplete(other.indStatscomplete());
  
  // keeping all from bim, but specifying N0etc need update
  newMat.setSnpStats(other.getSNPStats());
  newMat.setsnpStatscomplete(false);

  newMat.setMode(other.getMode());
}

template <typename SNPvectorClass, typename intVec>
SNPmatrix<SNPvectorClass> extractIndsfromSNPmatrixDisk(const SNPmatrix<SNPvectorClass> &other, const intVec &keep, std::string path_str) {
  SNPmatrix<SNPvectorClass> M;
  extractIndsfromSNPmatrixDisk(other, keep, path_str, M);
  return M;
}


#endif 
