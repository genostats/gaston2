#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr

#ifndef _datatype_
#define _datatype_

enum datatype { NONE, INT, DOUBLE, FLOAT, STRING };

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
