#include <string>
#include <stdexcept>
#include "isAutosome.h"

#ifndef _gaston_chrType_
#define _gaston_chrType_

enum chrType
{
  AUTOSOME,
  X,
  Y,
  MT,
  HAPLOTYPE,
  UNKNOWN
};

inline std::string chrTypeToString(chrType mo) {
  switch(mo) {
    case AUTOSOME:
      return std::string("AUTOSOME");
    case X:
      return std::string("X");
    case Y:
      return std::string("Y");
    case MT:
      return std::string("MT");
    case HAPLOTYPE:
      return std::string("HAPLOTYPE");
    case UNKNOWN:
      return std::string("UNKNOWN");
    default:
      throw std::runtime_error("something went horribly wrong (undefined chrType)");
  }
}

inline chrType stringTochrType(std::string str) {
  if(str == "AUTOSOME")  return chrType::AUTOSOME; 
  if(str == "X")         return chrType::X; 
  if(str == "Y")         return chrType::Y; 
  if(str == "MT")        return chrType::MT; 
  if(str == "HAPLOTYPE") return chrType::HAPLOTYPE; 
  if(str == "UNKNOWN")   return chrType::UNKNOWN; 
  throw std::runtime_error("something went horribly wrong (undefined chrType)");
}

inline chrType intToChrType(int chr) {
  if(isAutosome(chr))  return chrType::AUTOSOME;
  if(isX(chr))         return chrType::X;
  if(isY(chr))         return chrType::Y;
  if(isMt(chr))        return chrType::MT;
  return chrType::UNKNOWN;
}

#endif
