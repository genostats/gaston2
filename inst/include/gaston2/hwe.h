#include <cmath>
#include <limits>

#ifndef _gaston2_hwe_test_
#define _gaston2_hwe_test_

enum class HWEtest { chisquare, chisquare_yates, exact };

// workaround because partial specialization of templated function is not allowed
template<HWEtest test, typename scalar_t = double>
class hwe_f;

/* compute chisquare test for HWE, with or without yates continuity correction */

template<typename scalar_t, bool yates = true>
inline scalar_t hwe_chi_square(unsigned int a0, unsigned int a1, unsigned int a2) {
  unsigned int n = a0 + a1 + a2;
  if(n == 0) // no data
    return std::numeric_limits<scalar_t>::quiet_NaN();

  scalar_t p = (scalar_t) (2*a2+a1) / (scalar_t) (2*n);
  if(p == 0 || p == 1) return 1; // monomorphic SNP

  scalar_t e0 = (scalar_t) n * (1 - p) * (1 - p);
  scalar_t e1 = (scalar_t) n * 2 * p * (1 - p);
  scalar_t e2 = (scalar_t) n * p * p;

  scalar_t delta0 = (scalar_t) a0 - e0;
  scalar_t delta1 = (scalar_t) a1 - e1;
  scalar_t delta2 = (scalar_t) a2 - e2;

  if constexpr (yates) { // explicitely tell compiler that 'yates' is known at compile time
    // we don't use 0.5 as a correction but the min of 0.5 and abs(O - E) (as in R function for 2x2 chi² test)
    scalar_t cc = 0.5;
    delta0 = std::abs(delta0);
    delta1 = std::abs(delta1);
    delta2 = std::abs(delta2);

    cc = (cc < delta0)?cc:delta0;
    cc = (cc < delta1)?cc:delta1;
    cc = (cc < delta2)?cc:delta2;

    delta0 -= cc;
    delta1 -= cc;
    delta2 -= cc;
  }

  scalar_t chi2 = delta0*delta0/e0 + delta1*delta1/e1 + delta2*delta2/e2;

  // nice trick to compute the p value for chi² with 1 df
  return std::erfc(std::sqrt(chi2 / (scalar_t) 2.0));
}

template<typename scalar_t>
class hwe_f<HWEtest::chisquare, scalar_t> {
  public:
  inline scalar_t operator()(unsigned int a0, unsigned int a1, unsigned int a2) const {
    return hwe_chi_square<scalar_t, false>(a0, a1, a2);
  }
};

template<typename scalar_t>
class hwe_f<HWEtest::chisquare_yates, scalar_t> {
  public:
  inline scalar_t operator()(unsigned int a0, unsigned int a1, unsigned int a2) const {
    return hwe_chi_square<scalar_t, true>(a0, a1, a2);
  }
};


template<typename scalar_t>
class hwe_f<HWEtest::exact, scalar_t> {
  public:
  inline scalar_t operator()(unsigned int a0, unsigned int a1, unsigned int a2) const {
    unsigned int n = a0 + a1 + a2;
    if(n == 0)  // no data
      return std::numeric_limits<scalar_t>::quiet_NaN();

    if( (a0 == 0 && a1 == 0) || (a1 == 0 && a2 == 0)) // monomorphe
      return 1.0;

    if(a2 > a0) { // swap
      int tmp = a2;
      a2 = a0;
      a0 = tmp;
    }

    unsigned int m = 2*a2 + a1; // rare alleles
    unsigned int s = ((2*n - m)*m)/(2*n);
    if(s%2 != m%2) s++;	

    scalar_t grand_sum(1), small_sum(0);
    scalar_t target = 0;

    if(a1 < s) { // on calcule d'abord pour b1 < s
      // tail down and compute target
      int b2 = (m - s)/2;
      int b1 = s;
      int b0 = n - b1 - b2;
      scalar_t X = 1;
      bool over = false;
      while(b1 >= 2) {
        X *= (scalar_t) (b1*(b1 - 1)) / (scalar_t) (4*(1 + b0)*(1 + b2));
        grand_sum += X;
        if(b1 == a1 + 2) {
          target = X;
          over = true;
        }
        if(over) small_sum += X;
        b1 -= 2; b0++; b2++;
      }
      // tail up
      b2 = (m - s)/2;
      b1 = s;
      b0 = n - b1 - b2;
      X = 1;
      over = false;
      while(b1 <= m-2) {
        X *= (scalar_t) (4*(b0)*(b2)) / (scalar_t) ((b1 + 2)*(b1 + 1)) ;
        grand_sum += X;
        if(over)
          small_sum += X;
        else if(X <= target) {
          over = true;
          small_sum += X;
        }
        b1 += 2; b0--; b2--;
      }
    } else if(a1 > s) { // dans l'autre sens
      // tail up and compute target
      int b2 = (m - s)/2;
      int b1 = s;
      int b0 = n - b1 - b2;
      scalar_t X = 1;
      bool over = false;
      while(b1 <= m-2) {
        X *= (scalar_t) (4*(b0)*(b2)) / (scalar_t) ((b1 + 2)*(b1 + 1)) ;
        grand_sum += X;
        if(b1 + 2 == a1) {
          target = X;
          over = true;
        }
        if(over) small_sum += X;
        b1 += 2; b0--; b2--;
      }
      // tail down
      b2 = (m - s)/2;
      b1 = s;
      b0 = n - b1 - b2;
      X = 1;
      over = false;
      while(b1 >= 2) {
        X *= (scalar_t) (b1*(b1 - 1)) / (scalar_t) (4*(1 + b0)*(1 + b2));
        grand_sum += X;
        if(over)
          small_sum += X;
        else if(X <= target) {
          over = true;
          small_sum += X;
        }
        b1 -= 2; b0++; b2++;
      }
    } else { // a1 = s!! -> target = 1
      target = 1;
      // tail up
      int b2 = (m - s)/2;
      int b1 = s;
      int b0 = n - b1 - b2;
      scalar_t X = 1;
      bool over = false;
      while(b1 <= m-2) {
        X *= (scalar_t) (4*(b0)*(b2)) / (scalar_t) ((b1+2)*(b1+1)) ;
        grand_sum += X;
        if(over)
          small_sum += X;
        else if(X <= target) {
          over = true;
          small_sum += X;
        }
        b1 += 2; b0--; b2--;
      }
      // tail down
      b2 = (m - s)/2;
      b1 = s;
      b0 = n - b1 - b2;
      X = 1;
      over = false;
      while(b1 >= 2) {
        X *= (scalar_t) (b1*(b1 - 1)) / (scalar_t) (4*(1 + b0)*(1 + b2));
        grand_sum += X;
        if(over)
          small_sum += X;
        else if(X <= target) {
          over = true;
          small_sum += X;
        }
        b1 -= 2; b0++; b2++;
      }
    }
    if(target >= 1) small_sum += 1;
    return(small_sum/grand_sum);
  }
};
  
  
  
  
  
#endif
