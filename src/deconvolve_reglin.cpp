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

arma::mat decon_penalty_crossprod_cpp(int n_latent, const std::string& penalty) {
  if (penalty == "ridge") {
    return arma::eye<arma::mat>(n_latent, n_latent);
  }

  int differences = 2;
  if (penalty == "diff1") {
    differences = 1;
  } else if (penalty != "diff2") {
    stop("Unknown deconvolution penalty: " + penalty);
  }

  // Accumulate the short difference stencil directly. This is exactly D'D,
  // including its boundary entries, without forming D or multiplying dense
  // matrices that consist almost entirely of zeros.
  arma::mat crossprod(n_latent, n_latent, arma::fill::zeros);
  const double diff1[] = {-1.0, 1.0};
  const double diff2[] = {1.0, -2.0, 1.0};
  const double* stencil = differences == 1 ? diff1 : diff2;
  for (int row = 0; row < n_latent - differences; ++row) {
    for (int j = 0; j <= differences; ++j) {
      for (int k = 0; k <= differences; ++k) {
        crossprod(row + j, row + k) += stencil[j] * stencil[k];
      }
    }
  }
  return crossprod;
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

  arma::mat penalty_crossprod = decon_penalty_crossprod_cpp(n_latent, penalty);
  arma::mat best_activity(n_latent, n_signals, arma::fill::value(NA_REAL));

  bool tune_global = tune_by == "global";
  if (!tune_global && tune_by != "signal") {
    stop("tune_by must be 'global' or 'signal'");
  }
  if (gcv_rule != "min" && gcv_rule != "1se") {
    stop("gcv_rule must be 'min' or '1se'");
  }

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
    arma::mat conv_mat = decon_convolution_matrix_cpp(n_time, hrf);
    arma::mat fit_mat = prewhiten_gcv ? apply_ar1_whitening_cpp(conv_mat, ar1_rho) : conv_mat;
    arma::mat hth = fit_mat.t() * fit_mat;
    arma::mat hty = fit_mat.t() * y_fit;

    for (int li = 0; li < n_lambda; ++li) {
      double lambda = lambda_grid(li);
      arma::mat lhs = hth + lambda * penalty_crossprod;
      lhs.diag() += ridge_floor;

      arma::mat activity;
      double df = 0.0;
      if (!needs_scores) {
        best_activity = solve_sympd_fallback(lhs, hty);
        best_lambda_global = lambda;
        best_kernel_global = ki + 1;
        best_lambda_signal.fill(lambda);
        best_kernel_signal.fill(ki + 1);
        continue;
      }

      // Keep the band-aware solver for tuned/diagnostic fits. Benchmarks of
      // explicit dense factors and combined RHSs were slower for long runs.
      activity = solve_sympd_fallback(lhs, hty);
      df = arma::trace(solve_sympd_fallback(lhs, hth));
      arma::mat fitted = fit_mat * activity;
      arma::mat residual = y_fit - fitted;
      arma::rowvec rss_by_signal = arma::sum(arma::square(residual), 0);
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
          best_activity = activity;
          best_lambda_global = lambda;
          best_kernel_global = ki + 1;
        }
      } else {
        for (int si = 0; si < n_signals; ++si) {
          if (gcv_by_signal(si) < best_signal_score(si)) {
            best_signal_score(si) = gcv_by_signal(si);
            best_activity.col(si) = activity.col(si);
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
    arma::mat conv_mat = decon_convolution_matrix_cpp(n_time, hrf);
    arma::mat fit_mat = prewhiten_gcv ? apply_ar1_whitening_cpp(conv_mat, ar1_rho) : conv_mat;
    arma::mat hth = fit_mat.t() * fit_mat;
    arma::mat hty = fit_mat.t() * y_fit;
    arma::mat lhs = hth + selected_lambda * penalty_crossprod;
    lhs.diag() += ridge_floor;

    best_activity = solve_sympd_fallback(lhs, hty);
    best_lambda_global = selected_lambda;
    best_kernel_global = selected_kernel + 1;
  }

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
