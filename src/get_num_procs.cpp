#include <omp.h>

//[[Rcpp::export]]
int get_num_procs() {
  return omp_get_num_procs();
}

