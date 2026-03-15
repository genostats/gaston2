#include <stdexcept>
#include <numeric>
#include <cmath>
#include <limits>


// A function for the EM algorithm to compute LD 
// x must be an integer vector with size 9. 
template<typename scalar_t, typename intVec>
inline scalar_t LD_EM(intVec x, bool r_scale, scalar_t eps = std::sqrt(std::numeric_limits<scalar_t>::epsilon())) {

  if(x.size() != 9) throw std::runtime_error("In LD_EM, bad vector size");
 
  scalar_t deux_n = (scalar_t) (2 * std::accumulate(x.begin(), x.end(), 0));
  scalar_t fA = (scalar_t) (2*(x[0] + x[3] + x[6]) + x[1] + x[4] + x[7]) / deux_n;
  scalar_t fB = (scalar_t) (2*(x[0] + x[1] + x[2]) + x[3] + x[4] + x[5]) / deux_n;
  scalar_t fa = 1 - fA;
  scalar_t fb = 1 - fB;

  scalar_t beta = (scalar_t) x[4] / deux_n;
  scalar_t alpha0 = (scalar_t) (2*x[0] + x[1] + x[3]) / deux_n;
  scalar_t alpha1 = alpha0 + fa*fb - fA*fB;
  scalar_t alpha2 = fB - alpha0;
  scalar_t alpha3 = fA - alpha0;

  // this starting point works wonders when D' = 1 (corresponding to tau = 0 or tau = 1)
  // that is, the algorithm will converge in one interation
  // NOTE there might be several local maxima. We don't try to check this. This shouldn't
  // happen often when HW holds.
  //  scalar_t tau = (scalar_t) (2*x[0] + x[1] + x[3]) * (scalar_t) (2*x[8] + x[5] + x[7]);
  // tau /= tau + (scalar_t) (2*x[2] + x[1] + x[5]) * (scalar_t) (2*x[6] + x[3] + x[7]);
  scalar_t tau = 0.5; // The above trick gives bad results for some matrices when there's a large departure from HWE

  // some edge cases 
  if(tau != tau) { // NaN
    if(x[1] != 0 || x[8] != 0) 
      tau = 1.;
    else if(x[2] != 0 ||  x[6] != 0) 
      tau = 0.;
    else // whatever
      tau = 0.5; 
  }

  scalar_t tau_old = tau;
  do {
    scalar_t tau_beta = tau*beta;
    tau_old = tau;
    tau = (alpha0 + tau_beta)*(alpha1 + tau_beta);
    tau /= (tau + (alpha2 - tau_beta)*(alpha3 - tau_beta));
  } while(std::abs(tau-tau_old) > eps);

  scalar_t D = alpha0 + tau*beta - fA*fB;
  if(r_scale) 
    return D/std::sqrt(fa*fA*fb*fB);
  else 
    return D;
}

