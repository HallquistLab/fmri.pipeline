#include "fmri.pipeline.h"

namespace {

arma::mat decon_convolution_matrix_cpp(int n_time, const arma::vec& kernel) {
  int kernel_len = kernel.n_elem;
  int n_latent = n_time + kernel_len - 1;
  arma::mat conv_mat(n_time, n_latent, arma::fill::zeros);

  for (int ki = 0; ki < kernel_len; ++ki) {
    int lag = ki;
    for (int ti = 0; ti < n_time; ++ti) {
      int col = ti + kernel_len - 1 - lag;
      conv_mat(ti, col) = kernel(ki);
    }
  }

  return conv_mat;
}

arma::mat decon_penalty_band_cpp(int n_latent, int band, const std::string& penalty) {
  arma::mat p_band(band + 1, n_latent, arma::fill::zeros);
  double w0 = 0.0, w1 = 0.0, w2 = 0.0;
  if (penalty == "ridge") {
    w0 = 1.0;
  } else if (penalty == "diff1") {
    w1 = 1.0;
  } else if (penalty == "diff2") {
    w2 = 1.0;
  } else if (penalty == "matern1") {
    w0 = 0.10;
    w1 = 0.90;
  } else if (penalty == "sobolev") {
    w0 = 0.05;
    w1 = 0.45;
    w2 = 0.50;
  } else {
    stop("Unknown deconvolution penalty: " + penalty);
  }

  if (w0 > 0.0) {
    p_band.row(0) += w0;
  }
  if (w1 > 0.0) {
    const double diff1[] = {-1.0, 1.0};
    for (int row = 0; row < n_latent - 1; ++row) {
      for (int j = 0; j <= 1; ++j) {
        int col_j = row + j;
        for (int k = j; k <= 1; ++k) {
          int col_k = row + k;
          p_band(col_k - col_j, col_j) += w1 * diff1[j] * diff1[k];
        }
      }
    }
  }
  if (w2 > 0.0) {
    const double diff2[] = {1.0, -2.0, 1.0};
    for (int row = 0; row < n_latent - 2; ++row) {
      for (int j = 0; j <= 2; ++j) {
        int col_j = row + j;
        for (int k = j; k <= 2; ++k) {
          int col_k = row + k;
          p_band(col_k - col_j, col_j) += w2 * diff2[j] * diff2[k];
        }
      }
    }
  }
  return p_band;
}

void build_whitened_kernel_rows(
  const arma::vec& kernel,
  bool prewhiten,
  double rho,
  arma::vec& w0,
  arma::vec& w_rest,
  int& shift_rest
) {
  int kernel_len = kernel.n_elem;
  arma::vec rev_k(kernel_len);
  for (int i = 0; i < kernel_len; ++i) {
    rev_k(i) = kernel(kernel_len - 1 - i);
  }
  bool use_ar1 = prewhiten && std::isfinite(rho) && std::abs(rho) > std::numeric_limits<double>::epsilon();
  if (!use_ar1) {
    w0 = rev_k;
    w_rest = rev_k;
    shift_rest = 0;
  } else {
    w0 = std::sqrt(1.0 - rho * rho) * rev_k;
    w_rest.zeros(kernel_len + 1);
    w_rest.subvec(1, kernel_len) = rev_k;
    w_rest.subvec(0, kernel_len - 1) -= rho * rev_k;
    shift_rest = 1;
  }
}

arma::mat decon_hth_band_cpp(
  int n_time,
  int n_latent,
  int band,
  const arma::vec& w0,
  const arma::vec& w_rest,
  int shift_rest
) {
  arma::mat hth_band(band + 1, n_latent, arma::fill::zeros);
  int len0 = w0.n_elem;
  for (int j = 0; j < len0; ++j) {
    double wj = w0(j);
    for (int k = j; k < len0; ++k) {
      hth_band(k - j, j) += wj * w0(k);
    }
  }
  int len_r = w_rest.n_elem;
  arma::mat outer_r = w_rest * w_rest.t();
  for (int t = 1; t < n_time; ++t) {
    int base = t - shift_rest;
    for (int j = 0; j < len_r; ++j) {
      int col_j = base + j;
      for (int d = 0; d < len_r - j; ++d) {
        hth_band(d, col_j) += outer_r(j + d, j);
      }
    }
  }
  return hth_band;
}

// Compute hty_t = (fit_mat.t() * y_fit).t() of size n_signals x n_latent from y_t = y_fit.t()
arma::mat decon_hty_transposed_cpp(
  const arma::mat& y_t,
  int n_latent,
  const arma::vec& w0,
  const arma::vec& w_rest,
  int shift_rest
) {
  int n_signals = y_t.n_rows;
  int n_time = y_t.n_cols;
  arma::mat hty_t(n_signals, n_latent, arma::fill::zeros);

  int len0 = w0.n_elem;
  const arma::vec y0 = y_t.col(0);
  for (int j = 0; j < len0; ++j) {
    hty_t.col(j) += w0(j) * y0;
  }
  int len_r = w_rest.n_elem;
  for (int t = 1; t < n_time; ++t) {
    int base = t - shift_rest;
    const arma::vec yt = y_t.col(t);
    for (int j = 0; j < len_r; ++j) {
      hty_t.col(base + j) += w_rest(j) * yt;
    }
  }
  return hty_t;
}

// Compute residual rss_by_signal from activity_t (n_signals x n_latent) and y_t (n_signals x n_time)
arma::rowvec decon_rss_from_activity_t_cpp(
  const arma::mat& y_t,
  const arma::mat& activity_t,
  const arma::vec& w0,
  const arma::vec& w_rest,
  int shift_rest
) {
  int n_signals = y_t.n_rows;
  int n_time = y_t.n_cols;
  int len0 = w0.n_elem;
  int len_r = w_rest.n_elem;
  arma::vec rss(n_signals, arma::fill::zeros);
  arma::vec diff_t(n_signals);

  diff_t = y_t.col(0);
  for (int j = 0; j < len0; ++j) {
    diff_t -= w0(j) * activity_t.col(j);
  }
  rss += arma::square(diff_t);

  for (int t = 1; t < n_time; ++t) {
    int base = t - shift_rest;
    diff_t = y_t.col(t);
    for (int j = 0; j < len_r; ++j) {
      diff_t -= w_rest(j) * activity_t.col(base + j);
    }
    rss += arma::square(diff_t);
  }
  return rss.t();
}

// In-place lower-banded Cholesky factorization: A = L * L^T
// l_band has dimensions (band + 1) x n_latent, where l_band(d, j) = A(j + d, j)
bool chol_band_lower_inplace(arma::mat& l_band, int band) {
  int n = l_band.n_cols;
  for (int j = 0; j < n; ++j) {
    double diag_val = l_band(0, j);
    if (!std::isfinite(diag_val) || diag_val <= 0.0) {
      return false;
    }
    double l_jj = std::sqrt(diag_val);
    l_band(0, j) = l_jj;
    double inv_ljj = 1.0 / l_jj;
    int d_max = std::min(band, n - 1 - j);
    double* col_j = l_band.colptr(j);
    for (int d = 1; d <= d_max; ++d) {
      col_j[d] *= inv_ljj;
    }
    for (int k = 1; k <= d_max; ++k) {
      double l_jk = col_j[k];
      double* col_jk = l_band.colptr(j + k);
      for (int d = k; d <= d_max; ++d) {
        col_jk[d - k] -= col_j[d] * l_jk;
      }
    }
  }
  return true;
}

// Solve (L * L^T) * X = B where X_t = X^T and B_t = B^T have size (n_signals x n_latent)
void solve_band_chol_transposed(const arma::mat& l_band, int band, arma::mat& x_t) {
  int n = l_band.n_cols;
  int n_rhs = x_t.n_rows;
  // Forward substitution: L * Y = B
  for (int j = 0; j < n; ++j) {
    const double* l_col = l_band.colptr(j);
    double inv_ljj = 1.0 / l_col[0];
    double* x_j = x_t.colptr(j);
    for (int r = 0; r < n_rhs; ++r) {
      x_j[r] *= inv_ljj;
    }
    int d_max = std::min(band, n - 1 - j);
    for (int d = 1; d <= d_max; ++d) {
      double l_dj = l_col[d];
      double* x_jd = x_t.colptr(j + d);
      for (int r = 0; r < n_rhs; ++r) {
        x_jd[r] -= l_dj * x_j[r];
      }
    }
  }
  // Backward substitution: L^T * X = Y
  for (int j = n - 1; j >= 0; --j) {
    const double* l_col = l_band.colptr(j);
    double inv_ljj = 1.0 / l_col[0];
    double* x_j = x_t.colptr(j);
    int d_max = std::min(band, n - 1 - j);
    for (int d = 1; d <= d_max; ++d) {
      double l_dj = l_col[d];
      const double* x_jd = x_t.colptr(j + d);
      for (int r = 0; r < n_rhs; ++r) {
        x_j[r] -= l_dj * x_jd[r];
      }
    }
    for (int r = 0; r < n_rhs; ++r) {
      x_j[r] *= inv_ljj;
    }
  }
}

// Exact trace(A^{-1} * H^T H) using the Takahashi banded inverse recursion
double trace_band_takahashi(
  const arma::mat& l_band,
  const arma::mat& hth_band,
  int band,
  arma::mat& s_band
) {
  int n = l_band.n_cols;
  int stride = band + 1;
  double df = 0.0;

  for (int j = n - 1; j >= 0; --j) {
    const double* l_col = l_band.colptr(j);
    const double* h_col = hth_band.colptr(j);
    double* s_col = s_band.colptr(j);
    double inv_ljj = 1.0 / l_col[0];
    int d_max = std::min(band, n - 1 - j);

    for (int d = d_max; d >= 1; --d) {
      double sum = 0.0;
      const double* s_jd = s_band.colptr(j + d);
      // k < d: S(j + k, j + d) = s_band(d - k, j + k)
      for (int k = 1; k < d; ++k) {
        sum += l_col[k] * s_band.colptr(j + k)[d - k];
      }
      // k >= d: S(j + k, j + d) = s_band(k - d, j + d)
      for (int k = d; k <= d_max; ++k) {
        sum += l_col[k] * s_jd[k - d];
      }
      s_col[d] = -inv_ljj * sum;
    }

    double sum0 = 0.0;
    for (int k = 1; k <= d_max; ++k) {
      sum0 += l_col[k] * s_col[k];
    }
    s_col[0] = inv_ljj * (inv_ljj - sum0);

    double col_df = s_col[0] * h_col[0];
    for (int d = 1; d <= d_max; ++d) {
      col_df += 2.0 * s_col[d] * h_col[d];
    }
    df += col_df;
  }
  (void)stride;
  return df;
}

arma::mat solve_sympd_fallback(const arma::mat& lhs, const arma::mat& rhs) {
  arma::mat out;
  bool ok = arma::solve(out, lhs, rhs, arma::solve_opts::likely_sympd);
  if (!ok) {
    ok = arma::solve(out, lhs, rhs);
  }
  if (!ok) {
    stop("Regularized deconvolution linear solve failed");
  }
  return out;
}

double estimate_ar1_rho_cpp(const arma::mat& y) {
  int n_time = y.n_rows;
  int n_signals = y.n_cols;
  arma::vec rhos(n_signals, arma::fill::value(arma::datum::nan));

  for (int si = 0; si < n_signals; ++si) {
    arma::vec x0 = y.submat(0, si, n_time - 2, si);
    arma::vec x1 = y.submat(1, si, n_time - 1, si);
    x0 -= arma::mean(x0);
    x1 -= arma::mean(x1);
    double denom = std::sqrt(arma::accu(arma::square(x0)) * arma::accu(arma::square(x1)));
    if (denom > std::numeric_limits<double>::epsilon()) {
      rhos(si) = arma::dot(x0, x1) / denom;
    }
  }

  arma::uvec finite_idx = arma::find_finite(rhos);
  if (finite_idx.n_elem == 0) {
    return 0.0;
  }

  arma::vec finite_rhos = rhos.elem(finite_idx);
  double rho = arma::median(finite_rhos);
  if (rho > 0.95) rho = 0.95;
  if (rho < 0.0) rho = 0.0;
  return rho;
}

arma::mat apply_ar1_whitening_cpp(const arma::mat& x, double rho) {
  if (!std::isfinite(rho) || std::abs(rho) <= std::numeric_limits<double>::epsilon()) {
    return x;
  }

  arma::mat out = x;
  out.row(0) = std::sqrt(1.0 - rho * rho) * x.row(0);
  out.rows(1, x.n_rows - 1) = x.rows(1, x.n_rows - 1) - rho * x.rows(0, x.n_rows - 2);
  return out;
}

} // namespace

//' Rcpp backend for regularized linear deconvolution
//'
//' @name deconvolve_reglin_cpp
//' @param BOLDobs matrix of observed BOLD time series
//' @param kernels matrix of candidate HRF kernels, one per column
//' @param lambda_grid candidate regularization strengths
//' @param penalty regularization penalty
//' @param tune_by whether tuning is global or signal-specific
//' @param normalize whether to z-score columns of \code{BOLDobs}
//' @param demean whether to demean columns of \code{BOLDobs}
//' @param trim_kernel whether to trim the leading kernel support
//' @param ridge_floor small diagonal ridge for numerical stability
//' @param return_diagnostics whether to return diagnostics
//' @param prewhiten_gcv whether to use an AR(1)-prewhitened objective
//' @param prewhiten_rho optional AR(1) coefficient; \code{NA} estimates from data
//' @param gcv_rule GCV selection rule
//' @keywords internal
//'
//' @return deconvolved matrix or diagnostic list
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::RObject deconvolve_reglin_cpp(
  arma::mat BOLDobs,
  const arma::mat& kernels,
  const arma::vec& lambda_grid,
  const std::string& penalty,
  const std::string& tune_by,
  bool normalize,
  bool demean,
  bool trim_kernel,
  double ridge_floor,
  bool return_diagnostics,
  bool prewhiten_gcv,
  double prewhiten_rho,
  const std::string& gcv_rule
) {
  int n_time = BOLDobs.n_rows;
  int n_signals = BOLDobs.n_cols;
  int kernel_len = kernels.n_rows;
  int n_kernels = kernels.n_cols;
  int n_lambda = lambda_grid.n_elem;
  int n_latent = n_time + kernel_len - 1;

  if (n_time < 3) {
    stop("Need at least 3 timepoints for regularized linear deconvolution");
  }
  if (kernel_len < 2 || n_kernels < 1) {
    stop("At least one HRF kernel with length >= 2 is required");
  }
  if (n_lambda < 1) {
    stop("At least one lambda value is required");
  }

  arma::rowvec centers(n_signals, arma::fill::zeros);
  if (normalize || demean) centers = arma::mean(BOLDobs, 0);
  arma::rowvec scales(n_signals, arma::fill::ones);

  if (normalize) {
    scales = arma::stddev(BOLDobs, 0, 0);
    for (int i = 0; i < n_signals; ++i) {
      if (!std::isfinite(scales(i)) || scales(i) <= std::numeric_limits<double>::epsilon()) {
        scales(i) = 1.0;
      }
    }
    BOLDobs.each_row() -= centers;
    BOLDobs.each_row() /= scales;
  } else if (demean) {
    BOLDobs.each_row() -= centers;
  }

  double ar1_rho = 0.0;
  if (prewhiten_gcv) {
    ar1_rho = R_IsNA(prewhiten_rho) ? estimate_ar1_rho_cpp(BOLDobs) : prewhiten_rho;
    if (ar1_rho > 0.99 || ar1_rho < -0.99) {
      stop("prewhiten_rho must be between -0.99 and 0.99");
    }
  }
  arma::mat y_fit = prewhiten_gcv ? apply_ar1_whitening_cpp(BOLDobs, ar1_rho) : BOLDobs;
  arma::mat y_t = y_fit.t();

  bool tune_global = tune_by == "global";
  if (!tune_global && tune_by != "signal") {
    stop("tune_by must be 'global' or 'signal'");
  }
  if (gcv_rule != "min" && gcv_rule != "1se") {
    stop("gcv_rule must be 'min' or '1se'");
  }

  bool use_ar1 = prewhiten_gcv && std::isfinite(ar1_rho) &&
    std::abs(ar1_rho) > std::numeric_limits<double>::epsilon();
  int band_h = kernel_len - 1 + (use_ar1 ? 1 : 0);
  int band_p = (penalty == "ridge") ? 0 : ((penalty == "diff1" || penalty == "matern1") ? 1 : 2);
  int band = std::max(band_h, band_p);

  arma::mat penalty_band = decon_penalty_band_cpp(n_latent, band, penalty);
  arma::mat l_band(band + 1, n_latent);
  arma::mat s_band(band + 1, n_latent);
  arma::mat activity_t(n_signals, n_latent);
  arma::mat best_activity_t(n_signals, n_latent, arma::fill::value(NA_REAL));

  auto solve_dense_fallback = [&](const arma::vec& hrf, double lambda, arma::mat& out_activity_t, double& out_df, bool compute_df) {
    arma::mat conv_mat = decon_convolution_matrix_cpp(n_time, hrf);
    arma::mat fit_mat = prewhiten_gcv ? apply_ar1_whitening_cpp(conv_mat, ar1_rho) : conv_mat;
    arma::mat hth = fit_mat.t() * fit_mat;
    arma::mat hty = fit_mat.t() * y_fit;
    arma::mat lhs = hth;
    for (int j = 0; j < n_latent; ++j) {
      lhs(j, j) += lambda * penalty_band(0, j) + ridge_floor;
      int d_max = std::min(band, n_latent - 1 - j);
      for (int d = 1; d <= d_max; ++d) {
        double p_val = lambda * penalty_band(d, j);
        lhs(j + d, j) += p_val;
        lhs(j, j + d) += p_val;
      }
    }
    out_activity_t = solve_sympd_fallback(lhs, hty).t();
    if (compute_df) {
      out_df = arma::trace(solve_sympd_fallback(lhs, hth));
    }
  };

  double best_global_score = R_PosInf;
  arma::vec best_signal_score(n_signals, arma::fill::value(R_PosInf));
  arma::vec best_lambda_signal(n_signals, arma::fill::value(NA_REAL));
  arma::ivec best_kernel_signal(n_signals, arma::fill::value(NA_INTEGER));
  double best_lambda_global = NA_REAL;
  int best_kernel_global = NA_INTEGER;

  int n_tuning = n_kernels * n_lambda;
  arma::ivec tuning_kernel(n_tuning);
  arma::vec tuning_lambda(n_tuning);
  arma::vec tuning_df(n_tuning);
  arma::vec tuning_rss(n_tuning);
  arma::vec tuning_gcv(n_tuning);
  arma::vec tuning_gcv_se(n_tuning, arma::fill::zeros);
  int tuning_i = 0;
  // A single candidate needs no residuals, effective df, or GCV unless the
  // caller requests diagnostics. The inverse problem itself is unchanged.
  bool needs_scores = return_diagnostics || n_tuning > 1;

  for (int ki = 0; ki < n_kernels; ++ki) {
    arma::vec hrf = kernels.col(ki);
    arma::vec w0, w_rest;
    int shift_rest = 0;
    build_whitened_kernel_rows(hrf, prewhiten_gcv, ar1_rho, w0, w_rest, shift_rest);
    arma::mat hth_band = decon_hth_band_cpp(n_time, n_latent, band, w0, w_rest, shift_rest);
    arma::mat hty_t = decon_hty_transposed_cpp(y_t, n_latent, w0, w_rest, shift_rest);

    for (int li = 0; li < n_lambda; ++li) {
      double lambda = lambda_grid(li);
      l_band = hth_band + lambda * penalty_band;
      l_band.row(0) += ridge_floor;

      double df = 0.0;
      if (chol_band_lower_inplace(l_band, band)) {
        activity_t = hty_t;
        solve_band_chol_transposed(l_band, band, activity_t);
        if (!needs_scores) {
          best_activity_t = activity_t;
          best_lambda_global = lambda;
          best_kernel_global = ki + 1;
          best_lambda_signal.fill(lambda);
          best_kernel_signal.fill(ki + 1);
          continue;
        }
        df = trace_band_takahashi(l_band, hth_band, band, s_band);
      } else {
        solve_dense_fallback(hrf, lambda, activity_t, df, needs_scores);
        if (!needs_scores) {
          best_activity_t = activity_t;
          best_lambda_global = lambda;
          best_kernel_global = ki + 1;
          best_lambda_signal.fill(lambda);
          best_kernel_signal.fill(ki + 1);
          continue;
        }
      }

      arma::rowvec rss_by_signal = decon_rss_from_activity_t_cpp(y_t, activity_t, w0, w_rest, shift_rest);
      double denom = std::pow(n_time - df, 2.0);
      if (denom < std::numeric_limits<double>::epsilon()) {
        denom = std::numeric_limits<double>::epsilon();
      }
      arma::rowvec gcv_by_signal = rss_by_signal / denom;
      double mean_gcv = arma::mean(gcv_by_signal);
      double gcv_se = 0.0;
      if (n_signals > 1) {
        gcv_se = arma::stddev(gcv_by_signal) / std::sqrt(static_cast<double>(n_signals));
      }
      double total_rss = arma::accu(rss_by_signal);

      tuning_kernel(tuning_i) = ki + 1;
      tuning_lambda(tuning_i) = lambda;
      tuning_df(tuning_i) = df;
      tuning_rss(tuning_i) = total_rss;
      tuning_gcv(tuning_i) = mean_gcv;
      tuning_gcv_se(tuning_i) = gcv_se;
      ++tuning_i;

      if (tune_global) {
        if (mean_gcv < best_global_score) {
          best_global_score = mean_gcv;
          best_activity_t = activity_t;
          best_lambda_global = lambda;
          best_kernel_global = ki + 1;
        }
      } else {
        for (int si = 0; si < n_signals; ++si) {
          if (gcv_by_signal(si) < best_signal_score(si)) {
            best_signal_score(si) = gcv_by_signal(si);
            best_activity_t.row(si) = activity_t.row(si);
            best_lambda_signal(si) = lambda;
            best_kernel_signal(si) = ki + 1;
          }
        }
      }
    }
  }

  if (needs_scores && gcv_rule == "1se" && tune_global) {
    arma::uword min_idx = tuning_gcv.index_min();
    double threshold = tuning_gcv(min_idx) + tuning_gcv_se(min_idx);
    int selected_idx = min_idx;
    double selected_lambda = tuning_lambda(min_idx);
    double selected_gcv = tuning_gcv(min_idx);

    for (int ci = 0; ci < n_tuning; ++ci) {
      if (tuning_gcv(ci) <= threshold) {
        bool larger_lambda = tuning_lambda(ci) > selected_lambda;
        bool same_lambda_better = tuning_lambda(ci) == selected_lambda && tuning_gcv(ci) < selected_gcv;
        if (larger_lambda || same_lambda_better) {
          selected_idx = ci;
          selected_lambda = tuning_lambda(ci);
          selected_gcv = tuning_gcv(ci);
        }
      }
    }

    int selected_kernel = tuning_kernel(selected_idx) - 1;
    arma::vec hrf = kernels.col(selected_kernel);
    arma::vec w0, w_rest;
    int shift_rest = 0;
    build_whitened_kernel_rows(hrf, prewhiten_gcv, ar1_rho, w0, w_rest, shift_rest);
    arma::mat hth_band = decon_hth_band_cpp(n_time, n_latent, band, w0, w_rest, shift_rest);
    arma::mat hty_t = decon_hty_transposed_cpp(y_t, n_latent, w0, w_rest, shift_rest);
    l_band = hth_band + selected_lambda * penalty_band;
    l_band.row(0) += ridge_floor;
    if (chol_band_lower_inplace(l_band, band)) {
      best_activity_t = hty_t;
      solve_band_chol_transposed(l_band, band, best_activity_t);
    } else {
      double unused_df = 0.0;
      solve_dense_fallback(hrf, selected_lambda, best_activity_t, unused_df, false);
    }
    best_lambda_global = selected_lambda;
    best_kernel_global = selected_kernel + 1;
  }

  arma::mat best_activity = best_activity_t.t();
  if (trim_kernel) {
    best_activity = best_activity.rows(kernel_len - 1, best_activity.n_rows - 1);
  }

  Rcpp::NumericMatrix activity_out = Rcpp::wrap(best_activity);
  activity_out.attr("scaled:center") = Rcpp::wrap(centers);
  activity_out.attr("scaled:scale") = Rcpp::wrap(scales);

  if (!return_diagnostics) {
    return activity_out;
  }

  Rcpp::DataFrame tuning = Rcpp::DataFrame::create(
    Rcpp::Named("kernel") = tuning_kernel,
    Rcpp::Named("lambda") = tuning_lambda,
    Rcpp::Named("df") = tuning_df,
    Rcpp::Named("rss") = tuning_rss,
    Rcpp::Named("gcv") = tuning_gcv,
    Rcpp::Named("gcv_se") = tuning_gcv_se
  );

  if (tune_global) {
    return Rcpp::List::create(
      Rcpp::Named("activity") = activity_out,
      Rcpp::Named("lambda") = best_lambda_global,
      Rcpp::Named("kernel_index") = best_kernel_global,
      Rcpp::Named("tuning") = tuning,
      Rcpp::Named("penalty") = penalty,
      Rcpp::Named("tune_by") = tune_by,
      Rcpp::Named("gcv_rule") = gcv_rule,
      Rcpp::Named("prewhiten_gcv") = prewhiten_gcv,
      Rcpp::Named("ar1_rho") = ar1_rho,
      Rcpp::Named("scaled_center") = centers,
      Rcpp::Named("scaled_scale") = scales
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("activity") = activity_out,
    Rcpp::Named("lambda") = best_lambda_signal,
    Rcpp::Named("kernel_index") = best_kernel_signal,
    Rcpp::Named("tuning") = tuning,
    Rcpp::Named("penalty") = penalty,
    Rcpp::Named("tune_by") = tune_by,
    Rcpp::Named("gcv_rule") = gcv_rule,
    Rcpp::Named("prewhiten_gcv") = prewhiten_gcv,
    Rcpp::Named("ar1_rho") = ar1_rho,
    Rcpp::Named("scaled_center") = centers,
    Rcpp::Named("scaled_scale") = scales
  );
}
