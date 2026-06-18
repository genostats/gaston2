#include <stdexcept>
#include <string>

#include "isAutosome.h"

#ifndef _gaston_chrType_
#define _gaston_chrType_

enum class chrType {
  AUTOSOME,
  X,
  Y,
  MT,
  HAPLOTYPE,
  UNKNOWN     // if future evolutions, keep UNKNOWN as last element
              // it used in other places to know the number of values
};

inline std::string chrTypeToString(chrType mo) {
  switch (mo) {
    case chrType::AUTOSOME:
      return std::string("AUTOSOME");
    case chrType::X:
      return std::string("X");
    case chrType::Y:
      return std::string("Y");
    case chrType::MT:
      return std::string("MT");
    case chrType::HAPLOTYPE:
      return std::string("HAPLOTYPE");
    case chrType::UNKNOWN:
      return std::string("UNKNOWN");
    default:
      throw std::runtime_error("something went horribly wrong (undefined chrType)");
  }
}

inline chrType stringTochrType(std::string str) {
  if (str == "AUTOSOME") return chrType::AUTOSOME;
  if (str == "X") return chrType::X;
  if (str == "Y") return chrType::Y;
  if (str == "MT") return chrType::MT;
  if (str == "HAPLOTYPE") return chrType::HAPLOTYPE;
  if (str == "UNKNOWN") return chrType::UNKNOWN;
  throw std::runtime_error("something went horribly wrong (undefined chrType)");
}

inline chrType intToChrType(int chr) {
  if (isAutosome(chr)) return chrType::AUTOSOME;
  if (isX(chr)) return chrType::X;
  if (isY(chr)) return chrType::Y;
  if (isMt(chr)) return chrType::MT;
  return chrType::UNKNOWN;
}


// TODO vérifier l'utilité de cette fonction créée par Ju
//
// helper function to propagate the "right" chrType and
// also do a quick validity check
// called by the bindInds family
inline chrType wanted_chrType(chrType first, chrType second) {
  if (first != second) {
    if (first != chrType::UNKNOWN) {
      if (second != chrType::UNKNOWN) {
        throw std::logic_error("Error : chrTypes were incompatible");  // both are known butnot matching !
      }
      return first; // is defined and second is UNKNOWN
    } else {
      return second; // is defined but first isn't
    }
  } else { // both match
    return first;
  }
}

#endif
