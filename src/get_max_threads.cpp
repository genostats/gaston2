#include <omp.h>

//[[Rcpp::export]]
int get_max_threads() {
  return omp_get_max_threads();
}

