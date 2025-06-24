#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "SNPdosage.h"

#ifndef _snpdosagememory_
#define _snpdosagememory_

/**
 * @brief A class derived from SNPdosage abstract class
 *
 * if we want to replace with NA this works with R
 * data_.assign(nbInds_, std::nanf(""));
 * @emoji cry
 */
class SNPdosageMemory : public SNPdosage {
 private:
  std::vector<float> data_;

 public:
  // un constructeur, taille nbInds
  SNPdosageMemory(size_t nbInds) : SNPdosage(nbInds), data_(nbInds) {}

  // constructeur par copie
  SNPdosageMemory(const std::shared_ptr<SNPdosage> source) : SNPdosage(source->nbInds()),
                                                             data_(source->data(), source->data() + source->nbInds()) {}

  // constructeur par copie avec une sélection des individus,
  // par un vecteur d'index
  template <typename intVec>
  SNPdosageMemory(const std::shared_ptr<SNPdosage> source, intVec &keep) : SNPdosage(keep.size()) {
    const float *refdata = source->data();
    // check the keep (no idx to far or impossible) ?
    data_.assign(nbInds_, 0);  // creating of size == keep.size() with 0

    for (size_t i = 0; i < nbInds_; i++) {
      int ind_idx = keep[i];
      data_[i] = refdata[ind_idx];
    }
  }

  // constructeur concatenant 2 vecteurs
  SNPdosageMemory(const std::shared_ptr<SNPdosage> first, const std::shared_ptr<SNPdosage> second) : 
  SNPdosage(first->nbInds() + second->nbInds()), 
  data_(first->data(), first->data() + first->nbInds()) {
    data_.reserve(nbInds_);  // not to loose performance
    // now appending with second one
    std::copy(second->data(), second->data() + second->nbInds(), data_.data() + first->nbInds());
  }

  float *data() { return &data_[0]; }

  const float *data() const { return &data_[0]; }
};

#endif
