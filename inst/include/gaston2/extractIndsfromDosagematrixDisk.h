#include <cstring>
#include <iostream>
#include <memory>  // shared_ptr

#include "SNPmatrix.h"
#include "SNPdosageDisk.h"
#include "mio.hpp"

#ifndef _EXTRACTINDDOSMATRIX_DISK_
#define _EXTRACTINDDOSMATRIX_DISK_

// To be sure that I'm only extracting from a SNPmatrix with dosages,
// I removed the template of the SNPmatrix (previously SNPmatrix<SNPvectorClass>)

template <typename intVec>
void extractIndsfromDosagematrixDisk(const SNPmatrix<SNPdosage> &other, const intVec &keep, std::string path_str,
                                  SNPmatrix<SNPdosage> &newMat) {
  std::error_code error;

  /* FIRST : check if file exists, if it does, abort to not overwrite */
  const char * path_c = path_str.c_str();
  FILE *check = fopen(path_c, "rb");
  if (check) {
      fclose(check);
      throw std::runtime_error("The new matrix will overwrite an existing file, Aborting !");
  }

  /* THEN if it doesn't exist, go prepare it with the right size for mio */
  FILE *f = fopen(path_c, "wb");
  if (!f)
  {
    throw std::runtime_error("Failed to open file for writing");
  }

  const std::vector<std::shared_ptr<SNPdosage>> otherSNPs = other.getSNPs();
  
  int to_add = (keep.size() * sizeof(float)) * otherSNPs.size();
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
    newMat.push_back(std::make_shared<SNPdosageDisk<mio::access_mode::write>>(snp, file_ptr, SNPindex++, keep));
  }
  // keeping all SNPStats
  DataStruct og_snp_stats = other.getSNPStats();
  newMat.setSnpStats(og_snp_stats);
  //extract stats now and set stats_set_ to true
  DataStruct original_dt = other.getIndStats();
  newMat.setIndStats(DataStruct(original_dt, keep)); // for now does nothing else than fam file
  // but will when compute_IndStats will be implemented
}

template <typename intVec>
SNPmatrix<SNPdosage> extractIndsfromDosagematrixDisk(const SNPmatrix<SNPdosage> &other, const intVec &keep, std::string path_str) {
  SNPmatrix<SNPdosage> M;
  extractIndsfromDosagematrixDisk(other, keep, path_str, M);
  return M;
}

#endif
