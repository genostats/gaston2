#include "SNPvector.h"
#include "SNPvectorDisk.h"
#include <vector> // for omp reduction
#include <memory> // for shared_ptr
#include <stdexcept>  // for out of range exceptions

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
    if (nbSNPs() > 0 && nbInds() != v->nbInds())
    {
      std::cerr << "Pb loading SNP" << std::endl;
      throw std::out_of_range("Attempting to load a SNP with a different nb of individuals");
    }
    SNPs_.push_back(v);
    indStatsComputed_ = false; // si on ajoute des SNP les stats individuelles doivent être recalculées
  }

  /**
   * @brief get number of SNPs (see also nbSNPs)
   */
  unsigned int size() const { return SNPs_.size(); }
  /**
   * @brief get number of SNPs (see also size)
   */
  unsigned int nbSNPs() const { return SNPs_.size(); }
  /**
   * @brief get number of Inds
   */
  unsigned int nbInds() const {
    if(nbSNPs() == 0) {
      return 0;
    }
    return SNPs_[0]->nbInds();
  }

  // temporary func to test d°
  void deleteSNP()
  {
    SNPs_.pop_back();
  }

#pragma omp declare reduction(vec_int_plus : std::vector<int> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0))

  void compute_indStats()  {
    if(indStatsComputed_) return; // stats déjà calculées, on ne recalcule pas
                                  
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

        // check si byte entiere à lire
        int max_ind = 4;
        if (byte == nbBytes - 1)
        {
          int reste = nbInds % 4;
          if (reste)
            max_ind = reste; // if != 0, else stays at 4 to read full byte
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

    indStats_.push_back(N0s, "N0");
    indStats_.push_back(N1s, "N1");
    indStats_.push_back(N2s, "N2");
    indStats_.push_back(NAs, "NAs");

    indStatsComputed_ = true;
    // TODO : could add a checksum ?
  }

  // TODO à déplacer dans la classe DataStruct (fonction extractLines par exemple)
  DataStruct extract_indStats(const std::vector<size_t> &keep)
  {
    DataStruct filtered_stats;

    for (Column &col : indStats_.cols)
    {
      auto type = col.type();
      if (type == datatype::INT)
      {
        const auto &values = col.get<int>();
        std::vector<int> filtered;
        for (size_t idx : keep)
        {
          filtered.push_back(values->at(idx));
        }
        filtered_stats.push_back(filtered);
      }
      else if (type == datatype::FLOAT)
      {
        const auto &values = col.get<float>();
        std::vector<float> filtered;
        for (size_t idx : keep)
        {
          filtered.push_back(values->at(idx));
        }
        filtered_stats.push_back(filtered);
      }
      else if (type == datatype::DOUBLE)
      {
        const auto &values = col.get<double>();
        std::vector<double> filtered;
        for (size_t idx : keep)
        {
          filtered.push_back(values->at(idx));
        }
        filtered_stats.push_back(filtered);
      }
      else //find out how to add a custom type ? 
      {
        throw std::runtime_error("No an existing type yet");
      }
    }

    return filtered_stats;
  }

  // just a small helper f° to clrify extract_ind, not sure if usefull
  int read_ind(uint8_t byte, size_t ind_idx) {
    byte >>= (2 * (ind_idx % 4));
    return byte & 3; // extraction des bits correspondant
  }


  // creating a new SNPmatrix containing only the individuals with theiur index specified in "keep"
  // This new SNPmatrix can exists in memory (inMemory = true) or be filled with SNPvectorDisk. In the 
  // latter case, it's backed up by a file written with path and mapped with mio. 
  // In both cases the function returns the new matrix.
  SNPmatrix extract_ind(const std::vector<size_t> &keep, bool inMemory = true, std::string path = "default")
  {

    if (SNPs_.empty())
      throw std::runtime_error("Original matrix is empty !");

    if (inMemory && path != "default") throw std::runtime_error("Incorrect call, if the matrix is in memory, it doesn't need a path for a backup file");


    size_t ref_inds = SNPs_[0]->nbInds();
    for (int i : keep)
    {
      if (i < 0 || i >= ref_inds)
        throw std::out_of_range("No individuals with this index : " + std::to_string(i));
    }

    SNPmatrix new_matrix;
    size_t new_nb_inds = keep.size();

    if (inMemory)
    {
      for (const auto &snp : SNPs_)
      {
        // to feed to new matrix
        std::shared_ptr<SNPvectorMemory> new_snp = std::make_shared<SNPvectorMemory>(new_nb_inds);
        const uint8_t *refdata = snp->data();

        // Filling with 0, #https://www.cog-genomics.org/plink/1.9/formats#bed
        // If N is not divisible by four, the extra high-order bits in the last byte of each block are always zero.
        std::fill(new_snp->data(), new_snp->data() + new_snp->nbChars(), 0);

        for (size_t i = 0; i < new_nb_inds; ++i)
        {
          int ind_idx = keep[i];
          size_t currentChar = ind_idx / 4;    // index du byte
          int ind_gen = read_ind(refdata[currentChar], ind_idx);

          size_t new_byte = i / 4;
          size_t new_2bits = (i % 4) * 2;
          ind_gen <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte

          new_snp->data()[new_byte] |= ind_gen; // le || pour ne pas toucher au bits déjà set
        }
        new_matrix.SNPs_.push_back(new_snp);
      }
    }
    else
    {
      int nb_snps = SNPs_.size();
      std::error_code error;
      // check if file exists
      FILE *check = fopen(path.c_str(), "rb");
      if (check) {
          fclose(check);
          throw std::runtime_error("The new matrix will overwrite an existing file, Aborting !");
      }

      FILE *f = fopen(path.c_str(), "wb");
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
      int to_add = (new_nb_inds / 4 + ((new_nb_inds % 4 == 0u) ? 0 : 1)) * nb_snps + 3 ;
      if (fseek(f, to_add - 1, SEEK_SET) != 0)
      {
        fclose(f);
        throw std::runtime_error("Error when resizing file");
      }
      fputc(0, f); // this is what will size it up
      fclose(f);

      // will write using mio
      mio::mmap_sink file_ = mio::make_mmap_sink(path, 3, mio::map_entire_file, error);
      if (error)
      {
        std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
        throw std::runtime_error(errMsg);
      }
      char *data = file_.data(); // will be used to read
      int snp_idx = 0;

      for (const auto &snp : SNPs_)
      {
        // parcours de tous les snps de la matrices d'origine

        const uint8_t *refdata = snp->data();
        uint8_t newdata = 0; // exactly 1 byte

        for (size_t i = 0; i < new_nb_inds; ++i)
        {
          int ind_idx = keep[i];
          size_t currentChar = ind_idx / 4;    // index du byte
          int ind_gen = read_ind(refdata[currentChar], ind_idx);

          size_t new_byte = (i / 4);
          size_t new_2bits = (i % 4) * 2;

          ind_gen <<= (new_2bits); // shifter pour le mettre au bon endroit dans le nv byte
          newdata |= ind_gen;

          size_t snp_offset = snp_idx * ((new_nb_inds + 3) / 4); // ensure rounding up
          // once every 4 inds BYTE BY BYTE
          if (i % 4 == 3 || i == new_nb_inds - 1)
          {
            data[snp_offset + new_byte] = newdata;
            newdata = 0;
          }
        }
        // implies that I'm keeping at max the buffered SNP fully in memory before flushing
        file_.sync(error);
        if (error)
        {
          std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
          throw std::runtime_error(errMsg);
        }
        snp_idx++;
      }
      // this and the above call to `sync` have the same
      // effect as if the destructor had been invoked.
      file_.unmap(); // so I can open it again

      // no need to take into account offset here cos done in SNPVectorDisk constructor
      mio::mmap_source file_for_snps = mio::make_mmap_source(path, 0, mio::map_entire_file, error);
      std::shared_ptr<mio::mmap_source> file_ptr = std::make_shared<mio::mmap_source>(std::move(file_for_snps)); // don't know if good idea, creates a NEW sink
      for (size_t i = 0; i < nb_snps; i++)
      {
        std::shared_ptr<SNPVectorDisk<mio::access_mode::read>> snpVec(new SNPVectorDisk<mio::access_mode::read>(new_nb_inds, file_ptr, i)); // no mode array fr now
        new_matrix.push_back(snpVec);
      }
    }


    new_matrix.setIndStats(extract_indStats(keep)); //passing N0s N1s N2s and NAs for 
    return new_matrix;
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
    if (std::shared_ptr<SNPVectorDisk<mio::access_mode::read>> test = std::dynamic_pointer_cast<SNPVectorDisk<mio::access_mode::read>>(SNPs_[index]))
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

  const DataStruct getIndStats() const { return indStats_; }

  // TODO : see if by default possible ?
  // they need to be ordered !!!!
  void setIndStats(Column N0s, Column N1s, Column N2s, Column NAs)
  {
    indStats_.push_back(N0s);
    indStats_.push_back(N1s);
    indStats_.push_back(N2s);
    indStats_.push_back(NAs);
  }

  void setIndStats(DataStruct new_stats) {
    indStats_ = new_stats;
  }

private:
  DataStruct indStats_;
  std::vector<std::shared_ptr<SNPvector>> SNPs_;
  bool indStatsComputed_ = false;
};

#endif
