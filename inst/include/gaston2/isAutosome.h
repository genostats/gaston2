#include "gastonOptions.h"
#include <vector>

#ifndef _gaston_is_autosome_
#define _gaston_is_autosome_

template<typename intVector>
std::vector<bool> isAutosome(intVector & chr) {
  gastonOptions & opt = getGastonOptions();

  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) {
    test.push_back(opt.autosomes.find(c) != opt.autosomes.end());
  }
  return test;
}

template<typename intVector>
std::vector<bool> isX(intVector & chr) {
  gastonOptions & opt = getGastonOptions();

  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) {
    test.push_back(opt.x.find(c) != opt.x.end());
  }
  return test;
}

template<typename intVector>
std::vector<bool> isY(intVector & chr) {
  gastonOptions & opt = getGastonOptions();

  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) {
    test.push_back(opt.y.find(c) != opt.y.end());
  }
  return test;
}

template<typename intVector>
std::vector<bool> isMt(intVector & chr) {
  gastonOptions & opt = getGastonOptions();

  std::vector<bool> test;
  test.reserve(chr.size());
  for(int c : chr) {
    test.push_back(opt.mt.find(c) != opt.mt.end());
  }
  return test;
}

#endif
