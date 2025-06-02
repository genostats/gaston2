#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr

#ifndef _datatype_
#define _datatype_

enum datatype { NONE, INT, DOUBLE, FLOAT, STRING };

/***********************************************************/

inline std::string typeToString(datatype ty) {
  switch(ty) {
    case NONE:
      return std::string("NONE");
    case INT:
      return std::string("INT");
    case DOUBLE:
      return std::string("DOUBLE");
    case FLOAT:
      return std::string("FLOAT");
    case STRING:
      return std::string("STRING");
    default:
      throw std::runtime_error("unknown type");
  }
}

/***********************************************************/

template <typename T>
datatype whichType();

template<>
inline datatype whichType<int>() {
  return INT;
}

template<>
inline datatype whichType<double>() {
  return DOUBLE;
}

template<>
inline datatype whichType<float>() {
  return FLOAT;
}

template<>
inline datatype whichType<std::string>() {
  return STRING;
}

#endif
