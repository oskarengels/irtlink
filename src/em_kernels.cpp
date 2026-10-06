// src/em_kernels.cpp
//
// C++ kernels for the internal EM calibration engine

#include <Rcpp.h>
using namespace Rcpp;

// the E step lives in R/em-fit.R (em_estep)

// item M step, per-item Fisher scoring on the expected
// complete-data loglik (Hanson, 1998, Eq. 20). With
//   eta = a * theta - nu
// the score and the expected information are
//   score = ( sum_k (r - n P) theta_k , -sum_k (r - n P) )
//   info  = sum_k n P (1 - P) z z'   with   z = (theta_k, -1)
// The objective is concave, so scoring with step halving
// converges. est_a[j] = FALSE gives the scalar 1PL update.
// [[Rcpp::export]]
List em_mstep_items_cpp(NumericVector theta, NumericMatrix njk,
                        NumericMatrix rjk, NumericVector a,
                        NumericVector nu, LogicalVector est_a,
                        NumericVector lower_a, NumericVector upper_a,
                        NumericVector lower_nu, NumericVector upper_nu,
                        int inner_maxit, double inner_tol) {
  const int J = a.size(), K = theta.size();
  NumericVector a_new = clone(a), nu_new = clone(nu);
  for (int j = 0; j < J; ++j) {
    double aj = a_new[j], vj = nu_new[j];
    for (int it = 0; it < inner_maxit; ++it) {
      double s_a = 0, s_v = 0, i_aa = 0, i_av = 0, i_vv = 0;
      for (int k = 0; k < K; ++k) {
        const double n = njk(j, k);
        if (n <= 0) continue;
        const double r = rjk(j, k);
        double p = 1.0 / (1.0 + std::exp(-(aj * theta[k] - vj)));
        if (p < 1e-10) p = 1e-10;
        if (p > 1 - 1e-10) p = 1 - 1e-10;
        const double resid = r - n * p;
        const double w = n * p * (1 - p);
        s_a += resid * theta[k];
        s_v -= resid;
        i_aa += w * theta[k] * theta[k];
        i_av -= w * theta[k];
        i_vv += w;
      }
      double da = 0.0, dv = 0.0;
      if (est_a[j]) {
        const double det = i_aa * i_vv - i_av * i_av;
        if (det < 1e-12) break;
        da = ( i_vv * s_a - i_av * s_v) / det;
        dv = (-i_av * s_a + i_aa * s_v) / det;
      } else {
        if (i_vv < 1e-12) break;
        dv = s_v / i_vv;
      }
      double step = 1.0, aj1 = aj + da, vj1 = vj + dv;
      for (int h = 0; h < 12; ++h) {
        aj1 = aj + step * da;
        vj1 = vj + step * dv;
        if (aj1 >= lower_a[j] && aj1 <= upper_a[j] &&
            vj1 >= lower_nu[j] && vj1 <= upper_nu[j]) break;
        step *= 0.5;
      }
      aj1 = std::min(std::max(aj1, lower_a[j]), upper_a[j]);
      vj1 = std::min(std::max(vj1, lower_nu[j]), upper_nu[j]);
      const double delta = std::max(std::abs(aj1 - aj), std::abs(vj1 - vj));
      aj = aj1; vj = vj1;
      if (delta < inner_tol) break;
    }
    a_new[j] = aj; nu_new[j] = vj;
  }
  return List::create(_["a"] = a_new, _["nu"] = nu_new);
}
