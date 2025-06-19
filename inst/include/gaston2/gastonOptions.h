#include <vector>

#ifndef _gaston_options_
#define _gaston_options_

class gastonOptions {
  public:
  // defaut values = human plink conventions
  std::vector<int> autosomes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25};
  std::vector<int> x{23};
  std::vector<int> y{24};
  std::vector<int> mt{26};
};

// even if the function is declared inline,
// the static variable is shared between compilation units
// but don't declare it as static!
inline gastonOptions & getGastonOptions() {
  static gastonOptions opt;
  return opt;
}

inline void setGastonOptions(std::vector<int> & autosomes, std::vector<int> & x, std::vector<int> & y, std::vector<int> & mt) {
  gastonOptions & opt = getGastonOptions();
  opt.autosomes = autosomes;
  opt.x = x;
  opt.y = y;
  opt.mt = mt;
}

#endif
