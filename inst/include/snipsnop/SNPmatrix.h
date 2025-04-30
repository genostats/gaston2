#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector>
#include <memory>
#include <array>     // for compute stats function
#include <stdexcept> // for out of range exceptions
#include <filesystem> // for checking existence bedore writing
#include <fstream>      // std::ofstream


#include "Datastruct.h"

#include "SNPvectorMemory.h"

#include <omp.h>

#ifndef _snpmatrix_
#define _snpmatrix_

/**
 * @brief A class keeping in a vector pointers to SNPvectors
 * (either SNPvectorDisk or SNPvectorMemory)
 *
 */
class SNPmatrix
{
public:
  /**
   * @brief function to add a shared_ptr into the vector keeping the
   * SNPs of the matrix. /!\ will check if the SNP added is matching in size
   *
   */
  void push_back(std::shared_ptr<SNPvector> v)
  {
    // if at least one loaded, everySNPs must have same size
    if (SNPs_.size() && (SNPs_[0]->nbInds() != v->nbInds()))
    {
      std::cerr << "Pb loading SNP" << std::endl;
      throw std::out_of_range("Attempting to load a SNP with a different nb of individuals");
    }
    SNPs_.push_back(v);
  }

  int size() { return SNPs_.size(); }

  // temporary func to test d°
  void deleteSNP()
  {
    SNPs_.pop_back();
  }

  #pragma omp declare reduction(vec_int_plus : std::vector<int> : \
    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0))

void compute_stats()
  {
    int total = SNPs_.size();
    if (total == 0)
      throw std::out_of_range("No SNPs loaded into the SNPMatrix !");

    size_t nbInds = SNPs_[0]->nbInds();
    std::vector<int> unordered_stats(nbInds * 4, 0); // DO NOT USE nbCHars ! rounded up !

    #pragma omp parallel for reduction(vec_int_plus : unordered_stats)    
    for (int i = 0; i < total; i++)
    { // parcourt tous les SNPs
      auto snp = SNPs_[i];
      size_t nbBytes = snp->nbChars();

      // parcourt le SNP[i], byte by byte
      for (size_t byte = 0; byte < nbBytes; byte++)
      {
        uint8_t d = snp->data()[byte];
        size_t byteoffset = byte * 4;

        //check si byte entiere à lire
        int max_ind = 4;
        if (byte == nbBytes - 1)
        {
            int reste = nbInds % 4;
            if (reste) max_ind = reste; // if != 0, else stays at 4 to read full byte
        }

        for (int ind = 0; ind < max_ind; ind++)
        {
            unsigned int val_plink = snp->currentMode_[0][(d >> (2 * ind)) & 3];
            unordered_stats[(byteoffset + ind) * 4 + val_plink]++;
        }
      }
    }

    std::vector<int> vecN0s;
    std::vector<int> vecN1s;
    std::vector<int> vecN2s;
    std::vector<int> vecNAs;

    for (size_t ind = 0; ind < nbInds; ind++)
    {
      auto idxnzeros = ind * 4;
      vecN0s.push_back(unordered_stats[idxnzeros]);
      vecN1s.push_back(unordered_stats[idxnzeros + 1]);
      vecN2s.push_back(unordered_stats[idxnzeros + 2]);
      vecNAs.push_back(unordered_stats[idxnzeros + 3]);
    }

    Column N0s(vecN0s);
    Column N1s(vecN1s);
    Column N2s(vecN2s);
    Column NAs(vecNAs);

    stats_.push_back(N0s);
    stats_.push_back(N1s);
    stats_.push_back(N2s);
    stats_.push_back(NAs);
  }

  // WIP : NEED TO TEST THOROUGHLY, FIND A WAY FOR SNPVECTORDISK
  SNPmatrix extract_ind(const std::vector<size_t> &keep, bool inMemory= true)
  {

    if (SNPs_.empty())
      throw std::runtime_error("Original matrix is empty !");

    size_t ref_inds = SNPs_[0]->nbInds();
    for (int i : keep)
    {
      if (i < 0 || i >= ref_inds)
        throw std::out_of_range("No individuals with this index : " + std::to_string(i));
    }

    // TODO : displace later in loop 
    SNPmatrix new_matrix;
    size_t new_nb_inds = keep.size();

    if (inMemory) {
      for (const auto &snp : SNPs_)
      { // parcours de tous les snps de la matrices d'origine

        size_t i = 0;

          // to feed to new matrix
          std::shared_ptr<SNPvectorMemory> new_snp = std::make_shared<SNPvectorMemory>(new_nb_inds);
          const uint8_t *refdata = snp->data();

          // Filling with 0, #https://www.cog-genomics.org/plink/1.9/formats#bed
          // If N is not divisible by four, the extra high-order bits in the last byte of each block are always zero.
          std::fill(new_snp->data(), new_snp->data() + new_snp->nbChars(), 0);
          

          for (size_t i = 0; i < new_nb_inds; ++i)
          { // parcours de keep pour isoler lees individus
            int ind_idx = keep[i];
            // recherche de l'individus
            size_t currentChar = ind_idx / 4;    // index du byte
            uint8_t byte = refdata[currentChar]; // je récup la byte
            size_t current2bits = (ind_idx % 4); // where bits dans le byte
            byte >>= (2 * current2bits);

            int val1 = (byte & 3); // extraction des bits correspondant

            uint8_t val = val1; // pour avoir la bonne taille pour le ||

            size_t new_byte = i / 4;
            size_t new_2bits = (i % 4) * 2;

            val <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte

            // pê un masque possible ???
            new_snp->data()[new_byte] |= val; // le || pour ne pas toucher au bits déjà set
          }
          new_matrix.SNPs_.push_back(new_snp);
        } 
      return new_matrix;
    } else {
      std::cout << "You are trying to extract SNPs vectors on disk, which is currently being implemented\n";

      std::string path("extracted_disk_mat.txt");
      int nb_snps = SNPs_.size();
      std::error_code error;
      // TODO : see why filesystem not recognized
      // !! IMPORTANT to check if it exists without loading it !!
      //if (std::filesystem::exists(path)) throw std::runtime_error("The new matrix will overwrite an existing file, Aborting !");

      // TODO : check alternatives, will this load the full file in memory ? 
      std::ofstream file_test(path, std::ifstream::out); //will erase its content cos not opened in append
      if (file_test.bad()) throw std::runtime_error("This file does not exists\n");

      // Adding magic numbers to
      // Identify a bed file in SNP major mode
      try {
         file_test.put(108);
         file_test.put(27);
         file_test.put(1);
         int to_add =  (new_nb_inds/4 + ((new_nb_inds%4 == 0u)?0:1)) * nb_snps;
         // THIS WILL ALLOW TO RESIZE THE FILE, allocated to the full matrix future size
         file_test.seekp(to_add + 3, std::ios_base::end); // accounting for the 3 bytes already written
         file_test.write("", 1);// actually extend file
         file_test.flush();
         file_test.close();
      } catch (const std::out_of_range& e) {
        std::cerr << "Error when resizing file: " << e.what() << std::endl;
      }
      
      // maybe try to map to_add + 3 ?
    mio::mmap_sink file_ = mio::make_mmap_sink(path, 3, mio::map_entire_file, error);
    if (error) {
    std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
    throw std::runtime_error(errMsg); 
    }
    char* data = file_.data();// will be used to read

    for (const auto &snp : SNPs_)
    {
      // parcours de tous les snps de la matrices d'origine

      const uint8_t *refdata = snp->data();
      uint8_t newdata = 0; // a vector containing exactly 1 byte
      
        // WRITING BYTE BY BYTE (once every 4 individuals)
        for (size_t i = 0; i < new_nb_inds; ++i)
        { // parcours de keep pour isoler lees individus
          int ind_idx = keep[i];
          // recherche de l'individus
          size_t currentChar = ind_idx / 4;    // index du byte
          uint8_t byte = refdata[currentChar]; // je récup la byte
          size_t current2bits = (ind_idx % 4); // where bits dans le byte
          byte >>= (2 * current2bits);

          int val1 = (byte & 3); // extraction des bits correspondant
          uint8_t val = val1; // pour avoir la bonne taille pour le ||

          size_t new_byte = (i / 4); // no need to account for offset if given during 
          size_t new_2bits = (i % 4) * 2;

          val <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte

          newdata |= val;
          
          if (i%4 == 3 || i == new_nb_inds -1) { // once every 4 BYTE BY BYTE
            data[new_byte] = newdata; // put takes char and if it fails will throw an exception
            newdata = 0;
            //TODO : move sync one step lower
            file_.sync(error);
            if (error) { 
              std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
              throw std::runtime_error(errMsg);  
            }
          } 
        }
      }
      // this and the above call to `sync` have the same
      // effect as if the destructor had been invoked.
      file_.unmap(); // so I can open it again

      // no need to take into account offset here cos done in SNPVectorDisk constructor 
      mio::mmap_source file_for_snps = mio::make_mmap_source(path, 0, mio::map_entire_file, error);
      std::shared_ptr<mio::mmap_source> file_ptr = std::make_shared<mio::mmap_source>(std::move(file_for_snps)); // don't know if good idea, creates a NEW sink
      for(size_t i = 0; i < nb_snps; i++) {
        std::shared_ptr<SNPVectorDisk> snpVec(new SNPVectorDisk(new_nb_inds,file_ptr, i)); //no mode array fr now
        new_matrix.push_back(snpVec);
      }
      return new_matrix;
    }
  }

  /**
   * @brief A member function testing whether the SNP at index
   * is a SNPVectorDisk or a SNPVectorMemory
   *
   * @param index
   * @return true
   * @return false
   */
  bool onDisk(size_t index)
  { // TODO  check why was that written ????doesn't work, will see later
    if (std::shared_ptr<SNPVectorDisk> test = std::dynamic_pointer_cast<SNPVectorDisk>(SNPs_[index]))
      return true;
    return false; // test was a null ptr, not castable to a SNPVectorDisk
  }

  const std::vector<std::shared_ptr<SNPvector>> &getSNPs() const
  {
    return SNPs_;
  }

  const std::shared_ptr<SNPvector> &getSNP(size_t i) const
  {
    return SNPs_[i];
  }

  const DataStruct getStats() const { return stats_; }

private:
  DataStruct stats_;
  std::vector<std::shared_ptr<SNPvector>> SNPs_;
};

#endif