
// useful debugging macros
#ifndef SHOW
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << "\n";
#endif

#ifndef SHOWVEC
#define SHOWVEC(x) Rcpp::Rcout << #x << " = ";                \
  for(auto & _a_ : x) Rcpp::Rcout << _a_ << " ";                  \
  Rcpp::Rcout << "\n";
#endif
