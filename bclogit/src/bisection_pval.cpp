#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Compute HPD (Highest Posterior Density) interval from sorted samples.
// Returns a 2-element vector {lower, upper}.
static void hdi_sorted(const std::vector<double>& sorted_x, int n,
                        double level, double& lo, double& hi) {
  int window = static_cast<int>(std::floor(level * n));
  if (window < 1) window = 1;
  if (window >= n) {
    lo = sorted_x[0];
    hi = sorted_x[n - 1];
    return;
  }

  double min_width = sorted_x[window] - sorted_x[0];
  int best = 0;
  for (int i = 1; i + window < n; i++) {
    double w = sorted_x[i + window] - sorted_x[i];
    if (w < min_width) {
      min_width = w;
      best = i;
    }
  }
  lo = sorted_x[best];
  hi = sorted_x[best + window];
}

// Bisection p-value for a single coefficient's posterior samples.
// interval_type: 0 = HPD, 1 = CR (equal-tailed)
static double bisection_pval_single(const std::vector<double>& sorted_x,
                                     int n, int n_iter, int interval_type) {
  // Check if all samples are on one side of zero
  if (sorted_x[0] > 0.0 || sorted_x[n - 1] < 0.0) {
    return 0.0;
  }

  double low = 0.0, high = 1.0;
  for (int i = 0; i < n_iter; i++) {
    double mid = (low + high) / 2.0;
    if (mid < 0.001) mid = 0.001;
    if (mid > 0.999) mid = 0.999;

    double level = 1.0 - mid;
    bool in_interval;

    if (interval_type == 0) {
      // HPD
      double lo, hi;
      hdi_sorted(sorted_x, n, level, lo, hi);
      in_interval = (lo <= 0.0 && hi >= 0.0);
    } else {
      // Equal-tailed credible interval
      double tail_prob = mid / 2.0;
      // quantile using linear interpolation (Type 7, R default)
      double idx_lo = tail_prob * (n - 1);
      int il = static_cast<int>(std::floor(idx_lo));
      double fl = idx_lo - il;
      double q_lo = (il + 1 < n) ? sorted_x[il] * (1.0 - fl) + sorted_x[il + 1] * fl : sorted_x[il];

      double idx_hi = (1.0 - tail_prob) * (n - 1);
      int ih = static_cast<int>(std::floor(idx_hi));
      double fh = idx_hi - ih;
      double q_hi = (ih + 1 < n) ? sorted_x[ih] * (1.0 - fh) + sorted_x[ih + 1] * fh : sorted_x[ih];

      in_interval = (q_lo <= 0.0 && q_hi >= 0.0);
    }

    if (in_interval) {
      low = mid;
    } else {
      high = mid;
    }
  }
  return low;
}

// [[Rcpp::export(rng = false)]]
NumericVector calc_bisection_pvals_cpp(const NumericMatrix& beta_post,
                                        int n_iter,
                                        int interval_type) {
  int n_samps = beta_post.nrow();
  int n_coefs = beta_post.ncol();

  NumericVector pvals(n_coefs);
  std::vector<double> col_sorted(n_samps);

  for (int j = 0; j < n_coefs; j++) {
    for (int i = 0; i < n_samps; i++) {
      col_sorted[i] = beta_post(i, j);
    }
    std::sort(col_sorted.begin(), col_sorted.end());
    pvals[j] = bisection_pval_single(col_sorted, n_samps, n_iter, interval_type);
  }

  return pvals;
}
