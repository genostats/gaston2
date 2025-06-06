/**
 * @file readDosFileMemory.cpp
 * @author Hervé Perdry
 * @brief Main file with functions to interact with SNPvector and SNPmatrix
 * @date 2023-12-11
 * 
 */
 #ifndef READOSFILE_H
 #define READOSFILE_H

#include "SNPmatrix.h"
#include "SNPdosageMemory.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <memory>
#include "mio.hpp"
#include <cstring>
#include <Rcpp.h>

/**
 * @brief Reading a .dosf file, storing SNPs in a SNPmatrix returned
 * 
 * @param dosfile The name of the dos file
 * @param bimfile The name of the bim file
 * @param bamfile The name of the fam file
 * @param M reference to an originally empty SNPmatrix
 */

inline void readDosFileMemory(std::string dosfile, std::string bimfile, std::string famfile, SNPmatrix<SNPdosage> & M) {
  // first read bim and fam file, to determine size  
  M.readFamFile(famfile);
  M.readBimFile(bimfile);

  size_t n_snp = M.getSNPStats().nrow();
  size_t n_ind = M.getIndStats().nrow();

  // now ready to read dos file
  std::ifstream file(dosfile, std::ifstream::binary);
  if(!file.is_open()) {
    throw std::runtime_error("Cannot open file " + dosfile);
  }

  for(size_t i = 0; i < n_snp; i++) {
    //makes a shared_ptr on a vector of snips 
    std::shared_ptr<SNPdosageMemory> snpVec(new SNPdosageMemory(n_ind));
    float * data = snpVec->data();
    file.read(reinterpret_cast<char *>(data), n_ind * sizeof(float)); // check and be sure of what this does 
    M.push_back(snpVec);
  }
  file.close();
}

inline SNPmatrix<SNPdosage> readDosFileMemory(std::string bedfile, std::string bimfile, std::string famfile) {
  SNPmatrix<SNPdosage> M;
  readDosFileMemory(bedfile, bimfile, famfile, M);
  return M;
}

#endif /*READOSFILE_H*/
