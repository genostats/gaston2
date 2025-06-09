/**
 * @file readBedFileMemory.cpp
 * @author Hervé Perdry
 * @brief Main file with functions to interact with SNPvector and SNPmatrix
 * @date 2023-12-11
 * 
 */
#include "SNPmatrix.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <memory>
#include "mio.hpp"
#include <cstring>
#include <Rcpp.h>

#include "SNPdosage.h"
#include "readDosageFileMemory.h"
#include "readDosageFileDisk.h"

/**
 * @brief Reading a .dosf file, storing SNPs in a SNPmatrix returned
 * 
 * @param bedfile The name of the bed file
 * @param bimfile The name of the bim file
 * @param bamfile The name of the fam file
 * @param M reference to an originally empty SNPmatrix
 */

// R exported function
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosage>> readDosageFileMemory_(std::string bedfile, std::string bimfile, std::string famfile) {
  Rcpp::XPtr<SNPmatrix<SNPdosage>> pM(new SNPmatrix<SNPdosage>);
  readDosageFileMemory(bedfile, bimfile, famfile, *pM);
  return pM;
}

// R exported function
// [[Rcpp::export]]
Rcpp::XPtr<SNPmatrix<SNPdosage>> readDosageFileDisk_(std::string bedfile, std::string bimfile, std::string famfile) {
  Rcpp::XPtr<SNPmatrix<SNPdosage>> pM(new SNPmatrix<SNPdosage>);
  readDosageFileDisk(bedfile, bimfile, famfile, *pM);
  return pM;
}