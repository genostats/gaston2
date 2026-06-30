#include <omp.h>

//[[Rcpp::export]]
void set_num_threads(int num) {
  omp_set_num_threads(num);
}

