#include <set>

#ifndef _gaston_options_
#define _gaston_options_

class gastonOptions {
  public:
  // defaut values = human plink conventions
  std::set<int> autosomes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25};
  std::set<int> x{23};
  std::set<int> y{24};
  std::set<int> mt{26};
};

// a trick to have a static variable declared in a .h
//
// even if the function is declared inline, the linker will merge them 
// this way the static variable inside is shared between compilation units
inline gastonOptions & getGastonOptions() {
  static gastonOptions opt;
  return opt;
}

inline void setGastonOptions(std::set<int> & autosomes, std::set<int> & x, std::set<int> & y, std::set<int> & mt) {
  gastonOptions & opt = getGastonOptions();
  opt.autosomes = autosomes;
  opt.x = x;
  opt.y = y;
  opt.mt = mt;
}

#endif
