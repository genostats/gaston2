#include "SNPmatrix.h"
SNPmatrix readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp);
SNPmatrix readBedFileDisk(std::string filename, size_t n_ind, size_t n_snp, Mode mode = PLINK);
