#include <Rcpp.h>
#include "SNPmatrix.h"
#include "SNPvector.h"
#include "extractSNPmatrix.h"
#include "debug.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <omp.h>
//#include "Datastruct.h"

using namespace Rcpp;

SNPmatrix readBedFileMemory_old(std::string filename, size_t n_ind, size_t n_snp) {
  std::ifstream file(filename, std::ifstream::binary);
  if(!file.is_open()) {
    throw std::runtime_error("Cannot open file " + filename);
  }

  // check magic number
  char magic[3];
  file.read(magic, 3);
  if(magic[0] != 108 || magic[1] != 27) {
    throw std::runtime_error("Not a bed file");
  }
  if(magic[2] != 1) {
    throw std::runtime_error("Not a bed file in SNP major mode");
  } 

  SNPmatrix M;
  for(size_t i = 0; i < n_snp; i++) {
    //makes a shared_ptr on a vector of snips 
    std::shared_ptr<SNPvectorMemory> snpVec(new SNPvectorMemory(n_ind));
    size_t n = snpVec->nbChars();
    uint8_t * data = snpVec->data();
    file.read(reinterpret_cast<char *>(data), n);
    M.push_back(snpVec);
  }
  
  file.close();
  return M;
}


SNPmatrix readBedFileDisk_old(std::string path, size_t n_ind, size_t n_snp) {
  std::ifstream file_test(path, std::ifstream::binary);
  if (file_test.bad()) throw std::runtime_error("This file does not exists\n");
  std::error_code error;
  mio::mmap_source file_ = mio::make_mmap_source(path, 0, mio::map_entire_file, error);
  if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
  }
  std::shared_ptr<mio::mmap_source> file_ptr = std::make_shared<mio::mmap_source>(std::move(file_)); // don't know if good idea, creates a NEW sink
  // DON'T USE FILE_ FROM NOW ON COS IT'S NULL, WAS MOVED
  const char* data = reinterpret_cast<const char*>(file_ptr->data());// const necessary bcos read only 
  
  // check magic number
  char magic[3];
  for (int i = 0; i < 3; i++)
  {
    //std::cout << int(*data) << ' ';
    magic[i] = *(data + i);
  }
  //std::cout << '\n';
  if(magic[0] != 108 || magic[1] != 27) {
    throw std::runtime_error("Not a bed file");
  }
  if(magic[2] != 1) {
    throw std::runtime_error("Not a bed file in SNP major mode");
  }
  SNPmatrix M;

  auto file_offset = 3; // BCOS MAGIC BYTES 
  for(size_t i = 0; i < n_snp; i++) {
    std::shared_ptr<SNPvectorDisk<mio::access_mode::read>> snpVec(new SNPvectorDisk<mio::access_mode::read>(n_ind,file_ptr, i));
    //data should be taken by file_ptr

    size_t n = snpVec->nbChars(); // func inherited from SNPVec, gives back size used by SNP
    M.push_back(snpVec);
    // Increment the file offset based on the size of the data
    file_offset += n;
  }
  return M;
}


// [[Rcpp::export]]
IntegerVector test_readBedFileMemory(std::string filename, size_t n_ind, size_t n_snp, bool verbose = true)
{
  if (verbose) std::cout << " reading : " << filename << "\n n_ind : " << n_ind << "\n n_snp : " << n_snp << "\n";
  SNPmatrix M = readBedFileMemory_old(filename, n_ind, n_snp);
  std::vector<unsigned int> res;
  for (auto v : M.getSNPs())
  {
    res.push_back(v->sum());
  }
  return wrap(res);
}

// [[Rcpp::export]]
IntegerVector test_readBedFileDisk(std::string filename, size_t n_ind, size_t n_snp)
{
  SNPmatrix M = readBedFileDisk_old(filename, n_ind, n_snp);
  std::vector<unsigned int> res;
  for (auto v : M.getSNPs())
  {
    res.push_back(v->sum());
  }
  return wrap(res);
}

// [[Rcpp::export]]
IntegerVector test_delete(std::string filename, size_t n_ind, size_t n_snp)
{
  SNPmatrix M = readBedFileDisk_old(filename, n_ind, n_snp);
  std::vector<int> res;
  for (auto v : M.getSNPs())
  {
    res.push_back(v->sum());
  }
  SNPmatrix MM = readBedFileMemory_old(filename, n_ind, n_snp);
  for (auto v : MM.getSNPs())
  {
    M.push_back(v);
  }
  std::cout << "Is the 3rd SNP on disk ? (yes) : " << M.onDisk(3);
  std::cout << "\nAnd this one ? (no) : " << M.onDisk(n_snp + 3) << "\n";
  while (M.size() != 0)
  {
    M.deleteSNP();
  }
  std::cout << "Is the file destroyed before or after ?\n";
  return wrap(res);
}

IntegerVector loop_sum(SNPmatrix matrix)
{
  std::vector<int> res;
  for (auto v : matrix.getSNPs())
  {
    //== sum(503)
    v->compute_stats(); // to get stats full
    // v->compute_mu_sigma(); // needed to calculate mu and sigma
    res.push_back(v->sum());
  }
  return wrap(res);
}

/*
IntegerVector test_readModes(std::string filename, size_t n_ind, size_t n_snp)
{
  std::cout << " reading : " << filename << "\n n_ind : " << n_ind << "\n n_snp : " << n_snp << "\n";
  std::cout << " reading in Numeric Mode \n";
  SNPmatrix M = readBedFileMemory(filename, n_ind, n_snp, NUMERIC);
  IntegerVector res = loop_sum(M);
  std::cout << " reading in Centered (not implemented yet) Mode \n";
  M = readBedFileMemory(filename, n_ind, n_snp, CENTERED);
  IntegerVector res2 = loop_sum(M);
  for (auto i : res2)
    res.push_back(i);
  std::cout << " reading in Standardized (not implemented yet) Mode \n";
  M = readBedFileMemory(filename, n_ind, n_snp, STANDARDIZED_MU_SIGMA);
  IntegerVector res3 = loop_sum(M);
  for (auto i : res3)
    res.push_back(i);
  std::cout << " reading in PLINK Mode \n";
  M = readBedFileMemory(filename, n_ind, n_snp, PLINK);
  IntegerVector res4 = loop_sum(M);
  for (auto i : res4)
    res.push_back(i);
  return wrap(res);
}
*/

/************************
 *    Test SNP reading  *
 ************************/

// TODO : to test with numerous SNPs, so other rows in the bed matrix

const char file_hardcode[68] = "inst/extdata/LCT.bed";

// ON MEMORY, mean = 4.5 microseconds on average

// [[Rcpp::export]]
unsigned int test_performance_iterator_default(unsigned int n)
{
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503)
  { // cas où je dépasse SNP[0]
    nbSNPs = (n / 503) + 1;
    n = n % 503;
  } // parce que 503 individus dans le file hardcodé
  // std::cout << "Nb of SNPs: " << nbSNPs << " and nb of bits to read on last :" << n <<"\n";
  SNPmatrix M = readBedFileMemory_old(file_hardcode, (nbSNPs > 1) ? 503 : n, nbSNPs);
  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec = M.getSNP(i); // peut pas déréférencer là parce qu'instancie la classe abstraite SNPvector
    if (i == nbSNPs - 1)
    {
      S += vec->sum(n);
      // std::cout << "final S: " << S <<"\n";
      return S;
    }
    S += vec->sum();
    // std::cout << "S for now: " << S <<"\n";
  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_1(unsigned int n)
{
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503)
  { // cas où je dépasse SNP[0]
    nbSNPs = (n / 503) + 1;
    n = n % 503;
  }
  SNPmatrix M = readBedFileMemory_old(file_hardcode, (nbSNPs > 1) ? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec = M.getSNP(i);;
    if (i == nbSNPs - 1)
    {
      S += vec->sum(n);
      return S;
    }
    for (unsigned int a : *vec)
      S += a;
  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_2(unsigned int n)
{
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503)
  {
    nbSNPs = (n / 503) + 1;
    n = n % 503;
  }
  SNPmatrix M = readBedFileMemory_old(file_hardcode, (nbSNPs > 1) ? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec = M.getSNP(i);;
    if (i == nbSNPs - 1)
    {
      S += vec->sum(n);
      return S;
    }
    for (auto pa = vec->begin(); pa != vec->end(); ++pa)
      S += *pa;
  }
  return -1;
}

// ON DISK NOW, mean = 15 microseconds

// [[Rcpp::export]]
unsigned int test_performance_iterator_disk(unsigned int n)
{
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503)
  { // cas où je dépasse SNP[0]
    nbSNPs = (n / 503) + 1;
    n = n % 503;
  }
  SNPmatrix M = readBedFileDisk_old(file_hardcode, (nbSNPs > 1) ? 503 : n, nbSNPs);
  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec = M.getSNP(i);;
    if (i == nbSNPs - 1)
    {
      S += vec->sum(n);
      return S;
    }
    S += vec->sum();
  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_1d(unsigned int n)
{
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503)
  { // cas où je dépasse SNP[0]
    nbSNPs = (n / 503) + 1;
    n = n % 503;
  }
  SNPmatrix M = readBedFileMemory_old(file_hardcode, (nbSNPs > 1) ? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec =M.getSNP(i);;
    if (i == nbSNPs - 1)
    {
      S += vec->sum(n);
      return S;
    }
    for (unsigned int a : *vec)
      S += a;
  }
  return -1;
}

// [[Rcpp::export]]
unsigned int test_performance_iterator_2d(unsigned int n)
{
  int nbSNPs = 1;
  unsigned int S = 0;
  if (n > 503)
  {
    nbSNPs = (n / 503) + 1;
    n = n % 503;
  }
  SNPmatrix M = readBedFileDisk_old(file_hardcode, (nbSNPs > 1) ? 503 : n, nbSNPs);

  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec = M.getSNP(i);;
    if (i == nbSNPs - 1)
    {
      S += vec->sum(n);
      return S;
    }
    for (auto pa = vec->begin(); pa != vec->end(); ++pa)
      S += *pa;
  }
  return -1;
}

// TODO :

// [[Rcpp::export]]
IntegerVector test_snp_stats(unsigned int n)
{
  int nbSNPs = 1;
  if (n > 503)
  {
    nbSNPs = (n / 503) + 1;
    n = n % 503;
  }
  SNPmatrix M = readBedFileMemory_old(file_hardcode, (nbSNPs > 1) ? 503 : n, nbSNPs);

  std::vector<unsigned int> res(4, 0); // Initialize sum vector {0, 0, 0, 0}

  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec = M.getSNP(i);;
    vec->compute_stats();
    auto stats = vec->getStats();

    for (int j = 0; j < 4; j++)
    {
      res[j] += stats[j]; // Element-wise sum
    }
    // FOR now no way to stop on a part of the SNP
    if (i == nbSNPs - 1)
    {
      return wrap(res);
    }
  }
  return -1;
}

// [[Rcpp::export]]
IntegerMatrix test_snp_stats_all(int n_ind, int n_snp)
{

  SNPmatrix M = readBedFileMemory_old(file_hardcode, n_ind, n_snp);

  IntegerMatrix res(n_snp, 4);

  for (int i = 0; i < n_snp; i++)
  {
    auto vec = M.getSNP(i);
    vec->compute_stats();
    auto stats = vec->getStats();

    res(i, 0) = stats[0];
    res(i, 1) = stats[1];
    res(i, 2) = stats[2];
    res(i, 3) = stats[3];
  }
  return res; // retourne un vecteur avec les stats du SNP à la ligne correpondante
}

// [[Rcpp::export]]
IntegerVector test_snp_stats_d(unsigned int n)
{
  if (n > 503)
    n = 503;
  SNPmatrix M = readBedFileDisk_old(file_hardcode, n, 1);
  auto vec = M.getSNP(0);;
  vec->compute_stats();
  auto stats = vec->getStats();
  std::vector<unsigned int> res(stats, stats + 4);
  return wrap(res);
}

/*******************************
 *    Test SNP sum comparison  *
 ********************************/

// [[Rcpp::export]]
IntegerVector test_sums(unsigned int n)
{

  unsigned int sum_default = test_performance_iterator_default(n);
  auto sum_1 = test_performance_iterator_1(n);
  auto sum_2 = test_performance_iterator_2(n);
  auto sum_disk = test_performance_iterator_disk(n);
  auto sum_1d = test_performance_iterator_1d(n);
  auto sum_2d = test_performance_iterator_2d(n);

  if (std::set<unsigned int>{sum_default, sum_1, sum_2, sum_disk, sum_1d, sum_2d}.size() == 1)
  { // set rm duplicates, if size == 1 only 1 value
    std::cout << "All sums are equal!\n";
  }

  IntegerVector stats = test_snp_stats(n);
  return stats;
}

/*******************************
 *    Test mu and sigma  *
 ********************************/

// [[Rcpp::export]]
NumericMatrix test_mu_sigma(unsigned int n)
{
  int nbSNPs = n;

  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, nbSNPs);

  NumericMatrix res(2, nbSNPs);

  for (int i = 0; i < nbSNPs; i++)
  {
    auto vec = M.getSNP(i);
    // NEED TO COMPUTE BEFORE
    vec->compute_stats();
    // std::cout << "This is stats for SNP[" << i << "] : " << stats[0] << ", " << stats[1] << ", " << stats[2];

    vec->compute_mu_sigma();
    double sigma = vec->getSigma();
    double mu = vec->getMu();
    res(i, 0) = mu;
    res(i, 1) = sigma;
    // res.push_back(0);
    // std::cout << "These are stats for SNP[" << i << "], mu: " << mu << ", sigma: " << sigma << "\n";
  }

  return res;
}

// [[Rcpp::export]]
NumericVector test_modes_setsigma_one(int mode)
{
  std::vector<double> res;

  Mode mode_ = (Mode) mode;
  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, 607);
  M.setMode(mode_);
  for (auto v : M.getSNPs())
  {
    v->compute_stats();
    v->setSigma(1); // will call compute_mode also
    for (auto a : *v)
    {
      res.push_back(a);
    }
  }
  return wrap(res);
}

/********************************
 *           Test LD            *
 ********************************/

// [[Rcpp::export]]
NumericMatrix test_LD_square(int SNPnb1, int SNPnb2)
{
  if (SNPnb1 > 503 || SNPnb2 > 503)
  {
    std::cerr << "Please ju stop calling SNPs that don't exist \n"
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (SNPnb1 > SNPnb2)
  {
    std::cerr << "Please swap your SNPs to have a correct range \n"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  int nbSNPs = SNPnb2 - SNPnb1 + 1;
  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, 607); // should be good by loading aonly necessary snps
  NumericMatrix res(nbSNPs, nbSNPs);

  for (int i = 0; i < nbSNPs; i++)
  {
    SNPvector &snp1 = *M.getSNP(SNPnb1 + i);
    // first loading stats for mu and sigma:
    snp1.compute_stats();
    // snp1.compute_mu_sigma();

    for (int j = 0; j < nbSNPs; j++)
    {
      SNPvector &snp2 = *M.getSNP(SNPnb1 + j);
      if (j > i)
      { // aka where I haven't been before
        snp2.compute_stats();
        // snp2.compute_mu_sigma();
      }

       double r = snp1.LD(snp2); // parce que double par défaut dans template
       res(i, j) = r * r; // plus carré dans la fonction
    }
  }
  return res;
}

// [[Rcpp::export]]
SEXP test_contingency(int SNPnb1, int SNPnb2)
{
  if (SNPnb1 > 503 || SNPnb2 > 503)
  {
    std::cerr << "Please ju stop calling SNPs that don't exist \n"
              << std::endl;
    return 0;
  }
  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, 607); // should be good by loading aonly necessary snps ?

  IntegerVector res(9);
  SNPvector &snp1 = *M.getSNP(SNPnb1);
  SNPvector &snp2 = *M.getSNP(SNPnb2);
  snp1.contingency(snp2, res);
  res.attr("dim") = Dimension(3, 3);
  // IntegerMatrix m = as<IntegerMatrix>(res);
  // return m;
  return res;
}

// [[Rcpp::export]]
bool test_performance_stats_matrix()
{
  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, 607);
  M.compute_indStats();
  DataStruct st = M.getIndStats();

  //std::cout << "Stats loaded. Vector size: " << st.cols.size() << "\n";

  if (st.size() != 4) {
    std::cerr << "Error: expected 4 stat columns (N0, N1, N2, NA), got " << st.size() << "\n";
    return false;
  }

  std::ifstream file4("./inst/extdata/ind_counts.txt");
  int value;
  int i = 0;
  int j = 0;
  if (!file4)
  { 
    std::cout << "Pb with file\n";
    return false;
  }

  bool success = 1;

  // auto test1col = st.at(0);
  // std::cout << test1col.type() << "\n";
  // std::cout << test1col.get<int>() << "\n"; // THIS is a PTR to a vec of type int
  // std::cout << test1col.get<int>()->size();
  // std::cout << test1col.get<int>()->at(0);

  // sleep(2);

  while (file4 >> value)
  {
    //std::cout << " j = " << j << " and i = "<< i << "\n";
    //std::cout << "value " << value << "\n"; 
    if (i >= 503) {
      std::cerr << "Error: more inds in reference file than calculated !" << std::endl;
      return false;
    }
    
    auto computed = st.at(j).get<int>()->at(i);
    //std::cout << "this is computed val : " << computed << "\n";
    if ( computed != value)
    {
      std::cout << "Error: result_all_stats at (" << i << "," << j - 1 << ")  (line " << i + 1 << " and col n° " << j << " in the file) does not match reference value! " <<  st.at(j -1).get<int>()->at(i) << " != " << value<< std::endl;
      success = 0;
    }
    j++;

    //if (i > 500) sleep(1);

    if (j > 3) {
      i++;
      j = 0;
    }

  }
  if (!success) return false;
  return true;

}

// [[Rcpp::export]]
IntegerMatrix test_all_stats_matrix()
{
  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, 607); // should be good by loading aonly necessary snps
  //return wrap(M.compute_stats());

  M.compute_indStats();
  auto st = M.getIndStats();
  auto N0s = st.at(0).get<int>(); // a vec
  auto N1s = st.at(1).get<int>();
  auto N2s = st.at(2).get<int>();
  auto NAs = st.at(3).get<int>();

  IntegerMatrix res(503, 4);


  for (int i = 0; i < 503; i++)
  {
    int j = 0;
    res(i, j++) = N0s->at(i);
    res(i, j++) = N1s->at(i);
    res(i, j++) = N2s->at(i);
    res(i, j++) = NAs->at(i);
  }
  return res;
}

// //[[Rcpp::export]]
// void test_exception() {

// }

IntegerMatrix SNPmat_to_IntMat(SNPmatrix matrix)
{
  const std::vector<std::shared_ptr<SNPvector>> &SNPs = matrix.getSNPs();
  size_t n_snp = SNPs.size();
  if (n_snp == 0)
    return IntegerMatrix(0, 0);

  size_t n_ind = SNPs[0]->nbInds();
  IntegerMatrix mat(n_snp, n_ind);  // SNPs on rows, inds in cols

  for (size_t i = 0; i < n_snp; ++i) {
    auto snp = SNPs[i];
    size_t j = 0;
    for (auto it = snp->begin(); it != snp->end(); ++it) {
      mat(i, j++) = *it; // So read with PLINK mode !
    }
  }

  return mat;
}

// [[Rcpp::export]]
IntegerMatrix get_matrix_from_file(std::string path, size_t nbInds, size_t nbSNPs)
{
  SNPmatrix M = readBedFileDisk_old(path, nbInds, nbSNPs);
  IntegerMatrix x = SNPmat_to_IntMat(M);
  return x;
}

// [[Rcpp::export]]
void test_copyConstructor() 
{
      int nbInds = 503;
      int nbSNPs = 607;
      SNPmatrix M = readBedFileMemory_old(file_hardcode, nbInds, nbSNPs);

      FILE *f = fopen("/tmp/test.bed", "wb");
      if (!f)
      {
        throw std::runtime_error("Failed to open file for writing");
      }

      // Adding magic numbers to
      // Identify a bed file in SNP major mode
      fputc(108, f);
      fputc(27, f);
      fputc(1, f);
      // + 3 for the 3 magic bytes
      int to_add = (nbInds / 4 + ((nbInds % 4 == 0u) ? 0 : 1)) * nbSNPs + 3 ;
      if(fseek(f, to_add - 1, SEEK_SET) != 0)
      {
        fclose(f);
        throw std::runtime_error("Error when resizing file");
      }
      fputc(0, f); // this is what will size it up
      fclose(f);

      // will write using mio
      std::error_code error;
      mio::mmap_sink file_ = mio::make_mmap_sink("/tmp/test.bed", 0, mio::map_entire_file, error);
      if(error) throw std::runtime_error(error.message());

      std::shared_ptr<mio::mmap_sink> file_ptr = std::make_shared<mio::mmap_sink>(std::move(file_));
      // file_ can't be used anymore

      for(size_t snp = 0; snp < nbSNPs; snp++) {
        // j'ai pas besoin d'en garder trace pour ce test
        SNPvectorDisk<mio::access_mode::write>( M.getSNP(snp), file_ptr, snp);
      }
      std::cout << "created /tmp/test.bed, should be identical to " << file_hardcode << "\n";
}


//[[Rcpp::export]]
IntegerMatrix test_first_scnd_ind()
{
  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, 607);
  IntegerMatrix m(607, 3);
  int i = 0;
  for (auto v : M.getSNPs()) {
    auto it = v->begin();
    m(i, 0) = *(it);
    m(i, 1) = *(++it);
    m(i, 2) = *(++it); // read with mode
    i++;
  }
  return m;
}

//[[Rcpp::export]]
void set_num_thread(int num)
{
  omp_set_num_threads(num);
}








/****************************
 *        TESTSUITE         *
 ****************************/



#define GREEN "\033[1;32m"
#define RESET "\033[0m"
#define RED "\033[1;31m"

std::vector<std::string> tests_names = {
    "ReadBedFile from memory",
    "ReadBedFile from disk",
    "Computing stats",
    "Computing LD",
    "Computing values in centered mode",
    "Computing values in centered reduced mode",
    "Computing stats for all individuals", 
    "Computing extracted individuals", 
    "Extracting snps to a new matrix"
    // TODO : add a contingency test
  };

// helper to compare to doubles (surely very time consuming)
bool equal(double val1, double val2)
{
  double margin = 0.00000001;
  double res = val1 - val2;
  // std::cout << (res < margin) << " " << (-res < margin) << std::endl;

  return (res < margin) && (-res < margin);
}

// [[Rcpp::export]]
void testsuite(bool verbose = true)
{

  // if using valgrind, important
  set_num_thread(1);

  if (verbose) std::cout << "Using " << omp_get_max_threads() << " thread(s).\n";

  std::vector<int> total = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  // test_readBedFileMemory
  std::vector<int> expected = {804, 771, 982, 873, 399, 968, 976, 976, 976, 976, 976, 397, 873, 976, 873, 981, 976, 981, 843, 976, 804,
                               771, 464, 873, 843, 980, 399, 873, 398, 994, 771, 771, 843, 873, 771, 921, 981, 398, 795, 771, 995, 976,
                               771, 976, 921, 976, 803, 398, 398, 843, 968, 976, 398, 873, 991, 398, 976, 398, 756, 804, 368, 873, 398,
                               873, 771, 991, 873, 399, 873, 500, 976, 974, 976, 771, 974, 398, 398, 873, 803, 974, 771, 771, 771, 398,
                               771, 843, 994, 981, 976, 771, 873, 398, 398, 992, 873, 994, 769, 873, 974, 770, 397, 770, 770, 977, 397,
                               397, 991, 874, 977, 977, 770, 874, 875, 992, 769, 846, 977, 992, 977, 977, 397, 977, 847, 421, 976, 995,
                               977, 977, 769, 876, 396, 992, 847, 769, 769, 973, 977, 977, 977, 973, 976, 769, 976, 802, 860, 802, 973,
                               396, 876, 876, 876, 986, 979, 396, 845, 598, 598, 973, 973, 973, 875, 973, 363, 989, 977, 802, 925, 774,
                               887, 802, 993, 772, 801, 354, 772, 801, 801, 844, 974, 981, 363, 877, 831, 978, 988, 976, 994, 388, 976,
                               808, 389, 389, 965, 965, 857, 857, 857, 995, 857, 857, 806, 981, 995, 976, 390, 977, 976, 719, 719, 347,
                               975, 389, 389, 718, 992, 349, 718, 719, 994, 390, 965, 995, 872, 718, 719, 995, 983, 389, 791, 944, 872,
                               791, 945, 791, 965, 965, 872, 872, 976, 945, 976, 349, 791, 390, 976, 718, 945, 872, 390, 945, 389, 348,
                               873, 983, 389, 945, 995, 943, 965, 872, 793, 965, 977, 995, 392, 720, 347, 793, 994, 945, 719, 976, 719,
                               945, 945, 792, 975, 965, 945, 995, 792, 872, 389, 794, 792, 872, 872, 975, 720, 390, 720, 349, 718, 718,
                               720, 722, 965, 945, 995, 389, 793, 793, 347, 979, 976, 995, 975, 793, 346, 872, 793, 904, 995, 389, 995,
                               388, 792, 967, 991, 792, 792, 877, 793, 937, 377, 975, 874, 346, 388, 992, 874, 993, 990, 976, 926, 995,
                               792, 388, 791, 810, 792, 364, 784, 785, 845, 406, 407, 856, 974, 783, 753, 753, 753, 992, 855, 855, 907,
                               985, 992, 852, 977, 407, 855, 977, 993, 782, 855, 782, 977, 981, 855, 991, 964, 977, 977, 407, 977, 823,
                               977, 407, 993, 975, 989, 855, 855, 407, 407, 407, 855, 974, 974, 375, 782, 974, 782, 992, 407, 407, 782,
                               977, 782, 782, 976, 782, 973, 976, 976, 782, 856, 782, 407, 976, 782, 407, 990, 993, 783, 784, 406, 783,
                               856, 783, 403, 406, 856, 783, 856, 784, 980, 784, 855, 406, 406, 856, 985, 784, 976, 600, 976, 862, 976,
                               899, 793, 976, 631, 411, 411, 973, 977, 971, 973, 981, 970, 992, 987, 906, 413, 495, 973, 973, 411, 993,
                               419, 973, 419, 783, 985, 419, 994, 419, 419, 846, 494, 800, 419, 445, 783, 783, 970, 973, 970, 970, 623,
                               419, 783, 973, 973, 846, 783, 783, 419, 846, 783, 988, 783, 846, 419, 419, 419, 419, 995, 785, 846, 419,
                               781, 992, 781, 992, 419, 848, 848, 848, 781, 420, 848, 781, 419, 780, 848, 781, 419, 834, 418, 830, 793,
                               973, 849, 776, 993, 793, 973, 796, 833, 798, 798, 833, 798, 995, 798, 794, 977, 794, 975, 798, 832, 798,
                               794, 794, 797, 974, 417, 889, 838, 992, 798, 798, 797, 798, 784, 784, 780, 793, 990, 846, 784, 798, 417,
                               798, 832, 797, 801, 794, 784, 784, 974, 832, 780, 991, 417, 799, 832, 975, 798, 798, 417, 798, 794, 992,
                               795, 795, 784, 784, 386, 781, 798, 784, 798, 784, 417, 784, 795, 797, 800, 798, 976, 386, 994};

  IntegerVector result1 = test_readBedFileMemory("./inst/extdata/LCT.bed", 503, 607, verbose);
  IntegerVector result2 = test_readBedFileDisk("./inst/extdata/LCT.bed", 503, 607);
  if (result1.size() != expected.size())
    std::cerr << "Error: result from readBedFileMemory size does not match expected size!" << std::endl;
  else
    total[0] = 1;
  if (result2.size() != expected.size())
    std::cerr << "Error: result from readBedFileDisk size does not match expected size!" << std::endl;
  else
    total[1] = 1;
  for (size_t i = 0; i < expected.size(); i++)
  {
    if (result1[i] != expected[i])
    {
      std::cout << RED << "Error: result1 at index " << i << " does not match expected value!" << RESET << std::endl;
      total[0] = 0;
    }
    if (result2[i] != expected[i])
    {
      std::cout << RED << "Error: result2 at index " << i << " does not match expected value!" << RESET << std::endl;
      total[1] = 0;
    }
  }

  if (total[0])
  {
    if (verbose) std::cout << GREEN << "Tests for readBedFileMemory passed!" << RESET << std::endl;
  }
  if (total[1])
  {
    if (verbose) std::cout << GREEN << "Tests for readBedFileDisk passed!" << RESET << std::endl;
  }

  // snp_stats_all with ref file (snp_counts)
  std::vector<int> expected_stats;
  std::ifstream file("./inst/extdata/snp_counts.txt");
  int value;

  IntegerMatrix result_stats = test_snp_stats_all(503, 607);

  if (!file)
  {
    std::cerr << "Problem, failed to open reference file for test_snp_stats_all" << std::endl;
    exit(EXIT_FAILURE);
  }
  int i = 0;
  int j = 0;
  total[2] = 1;
  while (file >> value)
  {
    //std::cout << result_stats(i, j) << " and ref = " << value << std::endl;
    if (result_stats(i, j++) != value)
    {
      std::cout << RED << "Error: result_stats at (" << i << "," << j - 1 << ")  (line " << i + 1 << " and col n° " << j << " in the file) does not match reference value!" << RESET << std::endl;
      total[2] = 0;
    }
    if (i > 606)
      std::cerr << "Error: more snps in reference file than calculated !" << std::endl;
    if (j > 3)
    {
      i++;
      j = 0;
    }
    
  }
  if (total[2])
  {
    if (verbose) std::cout << GREEN << "Test for N0s N1s N2s and N3s on SNPs passed!" << RESET << std::endl;
  }

  // tests LD on all the snps
  std::ifstream file2("./inst/extdata/LD_ref.txt");
  double expected_LD;

  NumericMatrix result_LD = test_LD_square(0, 503);
  // if (result_LD == 0) goto conclusion;

  if (!file2)
  {
    std::cerr << "Problem, failed to open reference file for test_snp_stats_all" << std::endl;
    exit(EXIT_FAILURE);
  }
  i = 0;
  j = 0;
  total[3] = 1;

  while (file2 >> expected_LD)
  {
    //std::cout << result_LD(i, j) << " and ref = " << value2 << std::endl;
    if (!equal(result_LD(i, j), expected_LD))
    {
      std::cout << RED << "Error: result_LD at (" << i << "," << j << ")  (line " << i + 1 << " and col n° " << j << " in the file) does not match reference value!" << RESET << std::endl;
      total[3] = 0;
    }
    j++;
    if (i > 502)
    {
      std::cerr << "Error: more LD values in reference file than calculated !" << std::endl;
      total[3] = 0;
      break;
    }
    if (j > 502)
    {
      i++;
      j = 0;
    }
  }
  if (total[3])
  {
    if (verbose) std::cout << GREEN << "Test for LD values on all SNPs passed!" << RESET << std::endl;
  }

  // still need to check modes (how do i get them in gaston ??)
  NumericVector result_centered = test_modes_setsigma_one(1);
  NumericVector result_centered_reduced = test_modes_setsigma_one(2);
  std::ifstream file3("./inst/extdata/snps_centered.txt");
  // file generated from gaston : example("LCT"), x@sigma <- rep(1, ncol(x)), x@standardize_mu_sigma <- TRUE
  // !!! CHANGED EVERY NAs to 0, otherwise stopped the reading
  double expected_centered;

  if (!file3)
  {
    std::cerr << "Problem, failed to open reference file for test_snp_stats_all" << std::endl;
    exit(EXIT_FAILURE);
  }

  int max = result_centered.size();
  int cptr = 0;
  total[4] = 1;
  total[5] = 1;
  // std::cout << "This is max size  " << max << std::endl;
  while (file3 >> expected_centered)
  {

    if (cptr >= max)
    {
      // std::cout << "This is actual cptr  " << cptr << std::endl;
      std::cerr << "Error: more centered values in reference file than calculated !" << std::endl;
      total[4] = 0;
      break;
    }
    // std::cout << result_centered(cptr) << " and ref = " << expected_centered << std::endl;

    if (!equal(result_centered(cptr), expected_centered))
    {
      std::cout << RED << "Error: centered snp at line " << cptr + 1 << " in the file does not match computed value!" << RESET << std::endl;
      total[4] = 0;
    }
    if (!equal(result_centered_reduced(cptr), expected_centered))
    {
      std::cout << RED << "Error: centered reduced snp at line " << cptr + 1 << " in the file does not match computed value!" << RESET << std::endl;
      total[5] = 0;
    }
    cptr++;
  }
  if (total[4])
  {
    if (verbose) std::cout << GREEN << "Test for centered values on all SNPs passed!" << RESET << std::endl;
  }
  if (total[5])
  {
    if (verbose) std::cout << GREEN << "Test for centered_reduced values on all SNPs passed!" << RESET << std::endl;
  }

  // Now checking stats computed by individuals for full matrix
  // SNPmatrix::compute_stats with ref file (ind_counts)
  std::ifstream file4("./inst/extdata/ind_counts.txt");
  value = 0;

  if (!file4)
  {
    std::cerr << "Problem, failed to open reference file for test_stats_matrix" << std::endl;
    exit(EXIT_FAILURE);
  }

  i = 0;
  j = 0;
  total[6] = 1;
  IntegerMatrix result_stats_all = test_all_stats_matrix();

  while (file4 >> value)
  {
    if (i > 503)
      std::cerr << "Error: more inds in reference file than calculated !" << std::endl;
    
      //std::cout << result_stats_all(i * 4 +j) << " and ref = " << value << std::endl;

    if (result_stats_all(i,j++)!= value)
    {
      std::cout << RED << "Error: result_all_stats at (" << i << "," << j - 1 << ")  (line " << i + 1 << " and col n° " << j << " in the file) does not match reference value! " <<  result_stats_all(i,j-1) << " != " << value<< RESET << std::endl;
      total[6] = 0;
    }

    if (j > 3)
    {
      i++;
      j = 0;
    }

  }

  if (total[6])
  {
    if (verbose) std::cout << GREEN << "Test for N0s N1s N2s and N3s on inds passed!" << RESET << std::endl;
  }

  std::ifstream file5("./inst/extdata/extracted_snps.txt");
  value = 0;
  total[7] = 1;
  i = 0;
  j = 0;

  if (!file5)
  {
    std::cerr << "Problem, failed to open reference file extracted_snps.txt" << std::endl;
    exit(EXIT_FAILURE);
  }

  SNPmatrix M = readBedFileMemory_old(file_hardcode, 503, 607);
  M.compute_indStats();
  std::vector<size_t> to_keep = { 2, 6, 229, 230, 231, 232, 233, 234, 235, 236, 237 };
  SNPmatrix res_mat = extractIndsfromSNPmatrixMemory<std::vector<size_t>>(M, to_keep);
  IntegerMatrix res = SNPmat_to_IntMat(res_mat);
  IntegerMatrix res_disk = SNPmat_to_IntMat(extractIndsfromSNPmatrixDisk<std::vector<size_t>>(M, to_keep, "/tmp/extracted_mat_testsuite.bed"));
  
  IntegerMatrix from_file = get_matrix_from_file("/tmp/extracted_mat_testsuite.bed", 11, 607);

  while (file5 >> value)
  {
    if (i > 607) {
      std::cerr << "Error: more inds in reference file than calculated !" << std::endl;
      break;
    }
    
    if (from_file(i, j) != value) 
    {
      std::cout << RED << "Error: new extracted file with matrix at (" << i << "," << j << ") does not match reference value! " <<  from_file(i, j) << " != " << value << RESET << std::endl;
      total[7] = 0;
    }
       
    if (res(i,j)!= value)
    {
      std::cout << RED << "Error: new extracted matrix at (" << i << "," << j << ") does not match reference value! " <<  res(i,j) << " != " << value << RESET << std::endl;
      total[7] = 0;
    }
       
    if (res_disk(i,j++)!= value)
    {
      std::cout << RED << "Error: new extracted matrix on disk at (" << i << "," << j - 1 << ") does not match reference value! " <<  res(i,j - 1) << " != " << value << RESET << std::endl;
      total[7] = 0;
    }

    if ( i == 0 && j == 1) {
      DataStruct og_dt = M.getIndStats();
      Column og_N0 = og_dt.getColumn("N0");
      Column og_N1 = og_dt.getColumn("N1");
      Column og_N2 = og_dt.getColumn("N2");
      Column og_NA = og_dt.getColumn("NAs");

      DataStruct dt1 = res_mat.getIndStats();
      Column dt1_N0 = dt1.getColumn("N0");
      Column dt1_N1 = dt1.getColumn("N1");
      Column dt1_N2 = dt1.getColumn("N2");
      Column dt1_NA = dt1.getColumn("NAs");

      size_t t = 0;
      
      for (auto x : to_keep) {
        if ((og_N0.get<int>())[0][x] != (dt1_N0.get<int>())[0][t]) {
          std::cout << RED << "Error: stats for N0 on new extracted matrix (" << x << "th snp in the original one) are wrong !" << RESET << std::endl;
          total[8] = 0;
        }
        if ((og_N1.get<int>())[0][x] != (dt1_N1.get<int>())[0][t]) {
          std::cout << RED << "Error: stats for N1 on new extracted matrix (" << x << "th snp in the original one) are wrong !" << RESET << std::endl;
          total[8] = 0;
        }
        if ((og_N2.get<int>())[0][x] != (dt1_N2.get<int>())[0][t]) {
          std::cout << RED << "Error: stats for N2 on new extracted matrix (" << x << "th snp in the original one) are wrong !" << RESET << std::endl;
          total[8] = 0;
        }
        if ((og_NA.get<int>())[0][x] != (dt1_NA.get<int>())[0][t]) {
          std::cout << RED << "Error: stats for NA on new extracted matrix (" << x << "th snp in the original one) are wrong !" << RESET << std::endl;
          total[8] = 0;
        }
        t++;
      }
    }

    if (j > 10) // SNPs on rows, ind in cols
    {
      i++;
      j = 0;
    }
  }
  std::remove("/tmp/extracted_mat_testsuite.bed");
  if (total[7])
  {
    if (verbose) std::cout << GREEN << "Test for selected inds to extract in new matrix passed!" << RESET << std::endl;
  }

  std::vector<size_t> keep_idx = {0, 1, 2, 3, 100, 606};
  // now using the "interface" with C++, now sure how usefull it is,
  // allows to have it be more uniform I guess
  SNPmatrix extracted_snps = extractSNPsfromSNPmatrix(M, keep_idx);
  total[8] = 1;

  if (extracted_snps.getSNPs().size() != keep_idx.size()) {
    total[8] = 0;
    std::cout << RED << "Error: new extracted matrix for snps doesn't right nbr of new snps !" << RESET << std::endl;
  } else {
    for (size_t i = 0; i < keep_idx.size(); ++i) {
      if (extracted_snps.getSNPs()[i] != M.getSNPs()[keep_idx[i]]) {
        total[8] = 0;
        std::cout << RED << "Error: New SNPmatrix from extracted snps doesn't match at index " << i << "\n";
      }  
    }
  }
  if (total[8]) {
    if (verbose) std::cout << GREEN << "Test for selected snps to extract in new matrix passed!" << RESET << std::endl;
  }

conclusion:
  // std::count(vect.begin(), vect.end(), 0) should == 0 (all 1s) so !0 = true
  if (!std::count(total.begin(), total.end(), 0))
    std::cout << GREEN << total.size() << "/" << total.size() << " tests passed!" << RESET << std::endl;
  else
  {
    std::cout << RED << std::count(total.begin(), total.end(), 1) << "/" << total.size() << " tests passed..." << RESET << std::endl;
    for (int i = 0; i < total.size(); i++)
      if (!total[i])
        std::cout << RED << "Failed the " << tests_names[i] << " test..." << RESET << std::endl;
  }
}
