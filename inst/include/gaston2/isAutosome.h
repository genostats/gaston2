#include "gastonOptions.h"
#include <vector>

#ifndef _gaston_is_autosome_
#define _gaston_is_autosome_

// ------------- scalar functions ----------------------
inline bool isAutosome(int chr) {
  // static = shared with getGastonOptions which is called only once
  // -> can be called in loop without loss of performance
  static gastonOptions & opt = getGastonOptions();  
  return opt.autosomes.find(chr) != opt.autosomes.end();
}

inline bool isX(int chr) {
  static gastonOptions & opt = getGastonOptions();  
  return opt.autosomes.find(chr) != opt.x.end();
}

inline bool isY(int chr) {
  static gastonOptions & opt = getGastonOptions();  
  return opt.autosomes.find(chr) != opt.y.end();
}

inline bool isMt(int chr) {
  static gastonOptions & opt = getGastonOptions();  
  return opt.autosomes.find(chr) != opt.mt.end();
}

// ------------ vectorised functions -------------------
template<typename intVector>
std::vector<bool> isAutosome(intVector & chr) {
  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) test.push_back( isAutosome(c) );
  return test;
}

template<typename intVector>
std::vector<bool> isX(intVector & chr) {
  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) test.push_back( isX(c) );
  return test;
}

template<typename intVector>
std::vector<bool> isY(intVector & chr) {
  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) test.push_back( isY(c) );
  return test;
}

template<typename intVector>
std::vector<bool> isMt(intVector & chr) {
  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) test.push_back( isMt(c) );
  return test;
}

#endif
