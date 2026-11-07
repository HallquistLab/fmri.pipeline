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

// Marginal AR(1) innovations on the ORIGINAL acquisition clock. Missing rows
// are integrated out: a gap of d uses rho^d and its accumulated variance.
arma::mat hrf_whiten(const arma::mat& x, const arma::uvec& rows, double rho) {
  arma::mat out = x;
  if (rho == 0) return out;
  out.row(0) *= std::sqrt(1 - rho * rho);
  for (arma::uword i = 1; i < rows.n_elem; ++i) {
    double a = std::pow(rho, (double)(rows(i) - rows(i - 1)));
    out.row(i) = (x.row(i) - a * x.row(i - 1)) *
      std::sqrt((1 - rho * rho) / (1 - a * a));
  }
  return out;
}

arma::uvec hrf_rows(const Rcpp::IntegerVector& rows, arma::uword n) {
  if (rows.size() < 3) Rcpp::stop("At least three usable observations are required");
  arma::uvec out(rows.size());
  for (int i = 0; i < rows.size(); ++i) {
    if (rows[i] < 1 || rows[i] > (int)n || (i && rows[i] <= rows[i - 1]))
      Rcpp::stop("Observation rows must be increasing, unique, and within the run");
    out(i) = rows[i] - 1;
  }
  return out;
}

double hrf_f_tail(double extra, int rank, double rss, int df) {
  if (rank <= 0 || df < 2) return 1;
  if (rss <= 1e-14) return extra > 1e-12 ? 0 : 1;
  return R::pf(std::max(0.0, extra) * df / (rank * rss), rank, df, false, false);
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
    int batch_size = 2048, double rank_tol = 1e-8, double rho = 0,
    bool evidence = false, int fallback_index = 1) {
  if (bold.n_rows != nuisance.n_rows || bold.n_cols < 1 || rows.size() < 3 ||
      rows.size() > bold.n_rows || batch_size < 1 || designs.size() < 1 ||
      designs.size() > std::numeric_limits<int>::max() || rank_tol <= 0)
    Rcpp::stop("Invalid HRF scoring dimensions");
  if (!std::isfinite(rho) || std::abs(rho) >= 1 || fallback_index < 1 || fallback_index > designs.size())
    Rcpp::stop("Invalid AR coefficient or fallback index");
  arma::uvec keep = hrf_rows(rows, bold.n_rows);
  arma::mat z = nuisance.rows(keep);
  if (!z.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
  arma::mat qn = hrf_basis(hrf_whiten(z, keep, rho), rank_tol);
  int base_df = keep.n_elem - qn.n_cols;
  std::vector<arma::mat> qt(designs.size());
  Rcpp::IntegerVector ranks(designs.size());
  for (int k = 0; k < designs.size(); ++k) {
    arma::mat x = Rcpp::as<arma::mat>(designs[k]);
    if (x.n_rows != bold.n_rows || !x.is_finite()) Rcpp::stop("Invalid event design");
    x = hrf_whiten(x.rows(keep), keep, rho);
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
  int evidence_rows = evidence ? designs.size() : 0, evidence_cols = evidence ? bold.n_cols : 0;
  Rcpp::NumericMatrix task_p(evidence_rows, evidence_cols), fallback_p(evidence_rows, evidence_cols),
    incompatible_p(evidence_rows, evidence_cols);
  std::fill(task_p.begin(), task_p.end(), NA_REAL);
  std::fill(fallback_p.begin(), fallback_p.end(), NA_REAL);
  std::fill(incompatible_p.begin(), incompatible_p.end(), NA_REAL);
  struct Pair { int i, j, extra_rank, df; arma::mat q; };
  std::vector<Pair> pairs;
  if (comparable && evidence) for (int i = 0; i < designs.size(); ++i)
    for (int j = i + 1; j < designs.size(); ++j) {
      // Normalize only genuinely nonzero residual columns, so equal/subspace-
      // equivalent candidates never acquire rank from roundoff.
      arma::mat extra = qt[j] - qt[i] * (qt[i].t() * qt[j]);
      for (arma::uword c = 0; c < extra.n_cols; ++c)
        if (arma::norm(extra.col(c)) <= rank_tol) extra.col(c).zeros();
      arma::mat q = hrf_basis(extra, rank_tol);
      int d = q.n_cols;
      if (d && base_df - ranks[i] - d >= 2)
        pairs.push_back(Pair{i, j, d, base_df - ranks[i] - d, q});
    }
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
    y = hrf_whiten(y, keep, rho);
    if (qn.n_cols) y -= qn * (qn.t() * y);
    arma::rowvec sst = arma::sum(arma::square(y), 0);
    for (arma::uword j = 0; j < y.n_cols; ++j) {
      if (Rcpp::as<std::string>(status[start + j]) == "ok" && sst(j) <= rank_tol * rank_tol)
        status[start + j] = "nuisance_only";
      if (sst(j) > rank_tol * rank_tol) y.col(j) /= std::sqrt(sst(j));
    }
    arma::mat all_explained(designs.size(), y.n_cols);
    for (int k = 0; k < designs.size(); ++k) {
      arma::rowvec explained = arma::sum(arma::square(qt[k].t() * y), 0);
      for (arma::uword j = 0; j < y.n_cols; ++j) if (Rcpp::as<std::string>(status[start + j]) == "ok") {
        double r2 = std::max(0.0, std::min(1.0, explained(j)));
        score(k, start + j) = 1 - (1 - r2) * base_df / (base_df - ranks[k]);
        all_explained(k, j) = r2;
        if (evidence) {
          task_p(k, start + j) = hrf_f_tail(r2, ranks[k], 1 - r2, base_df - ranks[k]);
          fallback_p(k, start + j) = incompatible_p(k, start + j) = 1;
        }
      }
    }
    if (evidence) for (const Pair& pair : pairs) {
      arma::rowvec added = arma::sum(arma::square(pair.q.t() * y), 0);
      for (arma::uword j = 0; j < y.n_cols; ++j)
        if (Rcpp::as<std::string>(status[start + j]) == "ok") {
          double ei = all_explained(pair.i, j), ej = all_explained(pair.j, j);
          double both = std::min(1.0, std::max(std::max(ei, ej), ei + added(j)));
          double pi = hrf_f_tail(both - ei, pair.extra_rank, 1 - both, pair.df);
          double pj = hrf_f_tail(both - ej, pair.extra_rank, 1 - both, pair.df);
          // Only better-fitting competitors can exclude a candidate. Correct
          // for the whole competitor search later, across training runs too.
          if (ej > ei) incompatible_p(pair.i, start + j) =
            std::min(incompatible_p(pair.i, start + j), pi);
          if (ei > ej) incompatible_p(pair.j, start + j) =
            std::min(incompatible_p(pair.j, start + j), pj);
          if (pair.i == fallback_index - 1) fallback_p(pair.j, start + j) = pi;
          if (pair.j == fallback_index - 1) fallback_p(pair.i, start + j) = pj;
        }
    }
  }
  return Rcpp::List::create(Rcpp::Named("scores") = score,
    Rcpp::Named("signal_status") = status, Rcpp::Named("task_rank") = ranks,
    Rcpp::Named("nuisance_rank") = (int)qn.n_cols,
    Rcpp::Named("n_observations") = (int)keep.n_elem,
    Rcpp::Named("comparable") = comparable,
    Rcpp::Named("task_p") = task_p, Rcpp::Named("fallback_p") = fallback_p,
    Rcpp::Named("incompatible_p") = incompatible_p);
}

// REML on a fixed mean-model subspace, with a common factorization per rho and
// voxel blocks. The noise grid is deliberately small; no voxel optimizers.
// [[Rcpp::export]]
arma::mat hrf_noise_profile_cpp(const arma::mat& bold, const arma::mat& mean_design,
    const Rcpp::IntegerVector& rows, const arma::vec& rho_grid,
    int batch_size = 2048, double rank_tol = 1e-8) {
  arma::uvec keep = hrf_rows(rows, bold.n_rows);
  if (mean_design.n_rows != bold.n_rows || batch_size < 1) Rcpp::stop("Invalid noise design");
  arma::mat x = mean_design.rows(keep);
  if (!x.is_finite()) Rcpp::stop("Noise design must be finite on usable rows");
  arma::mat basis = hrf_basis(x, rank_tol);
  int df = keep.n_elem - basis.n_cols;
  arma::mat profile(rho_grid.n_elem, bold.n_cols, arma::fill::value(NA_REAL));
  if (df < 10) return profile;
  for (arma::uword k = 0; k < rho_grid.n_elem; ++k) {
    Rcpp::checkUserInterrupt();
    double rho = rho_grid(k);
    if (!std::isfinite(rho) || std::abs(rho) >= 1) Rcpp::stop("Invalid noise grid");
    arma::mat wx = hrf_whiten(basis, keep, rho), q, v;
    arma::vec s;
    if (!arma::svd_econ(q, s, v, wx)) Rcpp::stop("Noise design SVD failed");
    double logdet_design = 2 * arma::sum(arma::log(s));
    double logdet_w = .5 * std::log1p(-rho * rho);
    for (arma::uword i = 1; i < keep.n_elem; ++i) {
      double a = std::pow(rho, (double)(keep(i) - keep(i - 1)));
      logdet_w += .5 * (std::log1p(-rho * rho) - std::log1p(-a * a));
    }
    for (arma::uword start = 0; start < bold.n_cols; start += batch_size) {
      arma::uword end = std::min(bold.n_cols - 1, start + batch_size - 1);
      arma::mat y = bold.cols(start, end); y = y.rows(keep);
      arma::uvec valid(y.n_cols, arma::fill::ones);
      for (arma::uword j = 0; j < y.n_cols; ++j) {
        if (!y.col(j).is_finite()) { valid(j) = 0; y.col(j).zeros(); continue; }
        y.col(j) -= arma::mean(y.col(j));
        double norm = arma::norm(y.col(j));
        if (norm <= 0) valid(j) = 0;
        else y.col(j) /= norm;
      }
      y = hrf_whiten(y, keep, rho);
      arma::rowvec total = arma::sum(arma::square(y), 0), rss = total;
      if (q.n_cols) {
        arma::mat coefficients = q.t() * y;
        rss -= arma::sum(arma::square(coefficients), 0);
        // Avoid a second large BLAS multiplication for ordinary-SNR voxels.
        // Explicit residuals protect near-perfect fits from cancellation in
        // ||y||^2 - ||Q'y||^2 (including negative roundoff differences).
        for (arma::uword j = 0; j < y.n_cols; ++j)
          if (valid(j) && rss(j) < 1e-4 * total(j)) {
            arma::vec residual = y.col(j) - q * coefficients.col(j);
            rss(j) = arma::dot(residual, residual);
          }
      }
      for (arma::uword j = 0; j < y.n_cols; ++j) if (valid(j) && rss(j) > rank_tol * rank_tol)
        profile(k, start + j) = df * std::log(rss(j) / df) + logdet_design - 2 * logdet_w;
    }
  }
  return profile;
}

// Exact masked GLS inverse, jointly fitting unpenalized nuisance coefficients.
// Return only its action on measured rows; optional padded zeros are fixed.
// [[Rcpp::export]]
arma::mat hrf_inverse_operator_cpp(int n_time, const arma::vec& kernel,
    const Rcpp::IntegerVector& rows, double rho, const arma::mat& nuisance,
    double lambda, const std::string& penalty, double ridge_floor, int tail_padding = 0) {
  arma::uvec keep = hrf_rows(rows, n_time);
  if (kernel.n_elem < 2 || !kernel.is_finite() || nuisance.n_rows != (arma::uword)n_time ||
      !std::isfinite(rho) || std::abs(rho) >= 1 || lambda < 0 || ridge_floor < 0 || tail_padding < 0)
    Rcpp::stop("Invalid masked inverse parameters");
  arma::uword measured = keep.n_elem, history = kernel.n_elem - 1;
  if (tail_padding) keep = arma::join_cols(keep,
    arma::regspace<arma::uvec>(n_time, n_time + tail_padding - 1));
  arma::uword latent = n_time + tail_padding + history;
  arma::mat h(keep.n_elem, latent, arma::fill::zeros);
  for (arma::uword i = 0; i < keep.n_elem; ++i)
    for (arma::uword lag = 0; lag < kernel.n_elem; ++lag) h(i, keep(i) + history - lag) = kernel(lag);
  arma::mat w = hrf_whiten(arma::eye(keep.n_elem, keep.n_elem), keep, rho);
  h = w * h;
  arma::mat z(keep.n_elem, nuisance.n_cols, arma::fill::zeros);
  z.rows(0, measured - 1) = nuisance.rows(keep.head(measured));
  if (!z.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
  arma::mat q = hrf_basis(w * z, 1e-8);
  arma::mat rhs = w.cols(0, measured - 1);
  if (q.n_cols) { h -= q * (q.t() * h); rhs -= q * (q.t() * rhs); }
  arma::mat lhs = h.t() * h;
  if (penalty == "ridge") lhs.diag() += lambda;
  else {
    int order = penalty == "diff1" ? 1 : 2;
    if (penalty != "diff1" && penalty != "diff2") Rcpp::stop("Unknown inverse penalty");
    const double d1[] = {-1, 1}, d2[] = {1, -2, 1};
    const double* d = order == 1 ? d1 : d2;
    for (arma::uword i = 0; i + order < latent; ++i)
      for (int a = 0; a <= order; ++a) for (int b = 0; b <= order; ++b)
        lhs(i + a, i + b) += lambda * d[a] * d[b];
  }
  lhs.diag() += ridge_floor;
  arma::mat out;
  if (!arma::solve(out, lhs, h.t() * rhs, arma::solve_opts::likely_sympd))
    Rcpp::stop("Masked GLS inverse failed; increase ridge_floor");
  return out.rows(history, history + n_time - 1);
}
