#include "fmri.pipeline.h"

namespace {
arma::mat hrf_basis(const arma::mat& x, double tol) {
  if (x.n_cols == 0) return arma::mat(x.n_rows, 0);
  arma::mat scaled = x;
  for (arma::uword j = 0; j < scaled.n_cols; ++j) {
    double norm = arma::norm(scaled.col(j));
    if (norm > 0) scaled.col(j) /= norm;
  }
  arma::mat u, v;
  arma::vec s;
  if (!arma::svd_econ(u, s, v, scaled)) Rcpp::stop("HRF design SVD failed");
  if (s.n_elem == 0 || s(0) == 0) return arma::mat(x.n_rows, 0);
  arma::uword rank = arma::accu(s > tol * s(0));
  return rank ? u.cols(0, rank - 1) : arma::mat(x.n_rows, 0);
}
}

// Construct causal event designs on the acquisition grid. Fractional impulses
// use the same two-bin mass deposition as the impulse recovery diagnostic.
// Durations use fractional occupancy of [time_i, time_i + TR) bins.
// [[Rcpp::export]]
Rcpp::List hrf_event_design_cpp(int n_time, double tr, double offset,
    const arma::vec& onsets, const arma::vec& durations,
    const Rcpp::IntegerVector& conditions, int n_conditions,
    const arma::mat& kernels) {
  if (n_time < 3 || tr <= 0 || !std::isfinite(tr) ||
      onsets.n_elem != durations.n_elem || onsets.n_elem != conditions.size() ||
      n_conditions < 1 || kernels.n_rows < 2 || kernels.n_cols < 1)
    Rcpp::stop("Invalid HRF event design dimensions");
  // Include pre-run events whose causal response reaches observed samples.
  int history = kernels.n_rows - 1;
  arma::mat neural(n_time + history, n_conditions, arma::fill::zeros);
  for (arma::uword e = 0; e < onsets.n_elem; ++e) {
    if (!std::isfinite(onsets(e)) || !std::isfinite(durations(e)) || durations(e) < 0 ||
        conditions[e] < 1 || conditions[e] > n_conditions)
      Rcpp::stop("Invalid HRF event");
    double position = (onsets(e) - offset) / tr + history;
    int c = conditions[e] - 1;
    if (durations(e) == 0) {
      double lo = std::floor(position), frac = position - lo;
      if (lo >= 0 && lo < neural.n_rows) neural((int)lo, c) += 1 - frac;
      if (lo + 1 >= 0 && lo + 1 < neural.n_rows) neural((int)lo + 1, c) += frac;
    } else {
      double end = position + durations(e) / tr;
      int begin = (int)std::max(0.0, std::min((double)neural.n_rows, std::floor(position)));
      int stop = (int)std::max(0.0, std::min((double)neural.n_rows, std::ceil(end)));
      for (int t = begin; t < stop; ++t)
        neural(t, c) += std::max(0.0, std::min(end, t + 1.0) - std::max(position, (double)t));
    }
  }
  Rcpp::List designs(kernels.n_cols);
  for (arma::uword k = 0; k < kernels.n_cols; ++k) {
    arma::mat x(n_time, n_conditions, arma::fill::zeros);
    for (arma::uword lag = 0; lag < kernels.n_rows; ++lag)
      x += kernels(lag, k) * neural.rows(history - lag, history - lag + n_time - 1);
    designs[k] = x;
  }
  return designs;
}

// One SVD per candidate/run; block BLAS projections across all voxels. No
// voxelwise GLM objects, inverse matrices, or calls back into R are constructed.
// [[Rcpp::export]]
Rcpp::List hrf_score_run_cpp(const arma::mat& bold, const Rcpp::List& designs,
    const arma::mat& nuisance, const Rcpp::IntegerVector& rows,
    int batch_size = 2048, double rank_tol = 1e-8) {
  if (bold.n_rows != nuisance.n_rows || bold.n_cols < 1 || rows.size() < 3 ||
      rows.size() > bold.n_rows || batch_size < 1 || designs.size() < 1 ||
      designs.size() > std::numeric_limits<int>::max() || rank_tol <= 0)
    Rcpp::stop("Invalid HRF scoring dimensions");
  arma::uvec keep(static_cast<arma::uword>(rows.size()));
  for (int i = 0; i < rows.size(); ++i) {
    if (rows[i] < 1 || rows[i] > (int)bold.n_rows) Rcpp::stop("Invalid observation row");
    keep(i) = rows[i] - 1;
  }
  arma::mat z = nuisance.rows(keep);
  if (!z.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
  arma::mat qn = hrf_basis(z, rank_tol);
  int base_df = keep.n_elem - qn.n_cols;
  std::vector<arma::mat> qt(designs.size());
  Rcpp::IntegerVector ranks(designs.size());
  for (int k = 0; k < designs.size(); ++k) {
    arma::mat x = Rcpp::as<arma::mat>(designs[k]);
    if (x.n_rows != bold.n_rows || !x.is_finite()) Rcpp::stop("Invalid event design");
    x = x.rows(keep);
    arma::rowvec original_norm = arma::sqrt(arma::sum(arma::square(x), 0));
    if (qn.n_cols) x -= qn * (qn.t() * x);
    // A wholly absorbed task regressor must not become a spurious rank-one
    // model merely because relative SVD tolerance rescales roundoff to one.
    for (arma::uword j = 0; j < x.n_cols; ++j)
      if (arma::norm(x.col(j)) <= rank_tol * original_norm(j)) x.col(j).zeros();
    qt[k] = hrf_basis(x, rank_tol);
    ranks[k] = qt[k].n_cols;
  }
  bool comparable = base_df > 1 && ranks[0] > 0 && base_df > ranks[0] + 1;
  for (int k = 1; k < ranks.size(); ++k) comparable = comparable && ranks[k] == ranks[0];
  Rcpp::NumericMatrix score(static_cast<int>(designs.size()), bold.n_cols);
  std::fill(score.begin(), score.end(), NA_REAL);
  Rcpp::CharacterVector status(bold.n_cols, "unidentifiable_design");
  if (comparable) for (arma::uword start = 0; start < bold.n_cols; start += batch_size) {
    Rcpp::checkUserInterrupt();
    arma::uword end = std::min(bold.n_cols - 1, start + batch_size - 1);
    arma::mat y = bold.cols(start, end);
    y = y.rows(keep);
    for (arma::uword j = 0; j < y.n_cols; ++j) {
      if (!y.col(j).is_finite()) {
        status[start + j] = "nonfinite";
        y.col(j).zeros();
      } else {
        y.col(j) -= arma::mean(y.col(j));
        double norm = arma::norm(y.col(j));
        status[start + j] = norm > 0 ? "ok" : "constant";
        if (norm > 0) y.col(j) /= norm;
      }
    }
    if (qn.n_cols) y -= qn * (qn.t() * y);
    arma::rowvec sst = arma::sum(arma::square(y), 0);
    for (arma::uword j = 0; j < y.n_cols; ++j) {
      if (Rcpp::as<std::string>(status[start + j]) == "ok" && sst(j) <= rank_tol * rank_tol)
        status[start + j] = "nuisance_only";
      if (sst(j) > rank_tol * rank_tol) y.col(j) /= std::sqrt(sst(j));
    }
    for (int k = 0; k < designs.size(); ++k) {
      arma::rowvec explained = arma::sum(arma::square(qt[k].t() * y), 0);
      for (arma::uword j = 0; j < y.n_cols; ++j) if (Rcpp::as<std::string>(status[start + j]) == "ok") {
        double r2 = std::max(0.0, std::min(1.0, explained(j)));
        score(k, start + j) = 1 - (1 - r2) * base_df / (base_df - ranks[k]);
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("scores") = score,
    Rcpp::Named("signal_status") = status, Rcpp::Named("task_rank") = ranks,
    Rcpp::Named("nuisance_rank") = (int)qn.n_cols,
    Rcpp::Named("n_observations") = (int)keep.n_elem,
    Rcpp::Named("comparable") = comparable);
}
