#include <string>
#include <stdexcept>

#ifndef _gaston_mode_
#define _gaston_mode_

enum Mode
{
  RAW_VALUES,                   // .bed file : g = {0, 1, 3 and 2 = NA}
  CENTERED,              // g -= mu
  STANDARDIZED_MU_SIGMA, // g = (g-mu)/sd
  STANDARDIZED_P,        // g = (g - 2⁼p) / sqrt(2*p*(1-p)) with p = mu/2
  CUSTOM
};

inline std::string modeToString(Mode mo) {
  switch(mo) {
    case RAW_VALUES:
      return std::string("RAW_VALUES");
    case CENTERED:
      return std::string("CENTERED");
    case STANDARDIZED_MU_SIGMA:
      return std::string("STANDARDIZED_MU_SIGMA");
    case STANDARDIZED_P:
      return std::string("STANDARDIZED_P");
    case CUSTOM:
      return std::string("CUSTOM");
    default:
      throw std::runtime_error("unknown mode");
  }
}

inline Mode stringToMode(std::string str) {
  if(str == "RAW_VALUES") return Mode::RAW_VALUES; 
  if(str == "CENTERED") return Mode::CENTERED; 
  if(str == "STANDARDIZED_MU_SIGMA") return Mode::STANDARDIZED_MU_SIGMA; 
  if(str == "STANDARDIZED_P") return Mode::STANDARDIZED_P; 
  if(str == "CUSTOM") return Mode::CUSTOM; 
  throw std::runtime_error("unknown mode");
}

#endif
