#include "fmri.pipeline.h"

namespace {
arma::mat hrf_basis(const arma::mat& x, double tol) {
  if (x.n_cols == 0) return arma::mat(x.n_rows, 0);
  arma::mat scaled = x;
  for (arma::uword j = 0; j < scaled.n_cols; ++j) {
    double norm = arma::norm(scaled.col(j));
    if (norm > 0) scaled.col(j) /= norm;
  }
  if (scaled.n_rows > 2 * scaled.n_cols && scaled.n_cols >= 4) {
    arma::mat q_qr, r_qr, u_r, v_r;
    arma::vec s;
    if (arma::qr_econ(q_qr, r_qr, scaled) && arma::svd(u_r, s, v_r, r_qr)) {
      if (s.n_elem == 0 || s(0) == 0) return arma::mat(x.n_rows, 0);
      arma::uword rank = arma::accu(s > tol * s(0));
      return rank ? (q_qr * u_r.cols(0, rank - 1)) : arma::mat(x.n_rows, 0);
    }
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
void hrf_whiten_coefs(const arma::uvec& rows, double rho, arma::vec& d, arma::vec& s,
                      bool unit_marginal = false) {
  arma::uword m = rows.n_elem;
  d.ones(m);
  s.zeros(m);
  if (rho == 0 || m == 0) return;
  double s0 = unit_marginal ? 1.0 : std::sqrt(std::max(1.0 - rho * rho, 1e-8));
  d(0) = s0;
  for (arma::uword i = 1; i < m; ++i) {
    double a = std::pow(rho, (double)(rows(i) - rows(i - 1)));
    double mult = s0 / std::sqrt(std::max(1.0 - a * a, 1e-12));
    d(i) = mult;
    s(i) = -a * mult;
  }
}

void hrf_whiten_inplace(arma::mat& x, const arma::vec& d, const arma::vec& s, double rho) {
  if (rho == 0 || x.n_cols == 0 || x.n_rows == 0) return;
  const arma::uword m = x.n_rows;
  const double* dp = d.memptr();
  const double* sp = s.memptr();
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    double* col = x.colptr(j);
    for (arma::uword i = m - 1; i >= 1; --i) {
      col[i] = dp[i] * col[i] + sp[i] * col[i - 1];
    }
    col[0] *= dp[0];
  }
}

// Unified AR(1) / cascaded AR(p) whitener on observed acquisition indices `keep`.
// Matches `hrf_whitening_matrix(rows, ar_vec, unit_marginal)` to machine precision
// in O(m * p) time and O(p) memory on BOTH contiguous and motion-scrubbed grids:
// - Pure AR(1): O(m) bidiagonal pass.
// - Interior AR(p) rows (where p preceding samples on the clock are observed):
//   exact p-tap causal innovation filter.
// - Run-head & post-gap AR(p) boundary rows: exact local Cholesky back to the
//   nearest preceding p-sample Markov barrier (d-separating all earlier samples).
struct ArWhitener {
  struct BoundaryRow {
    arma::uword b0 = 0;
    arma::vec w_unit;
  };

  bool is_ar_p = false;
  bool use_dense_fallback = false;
  bool unit_marginal = false;
  double rho = 0.0;
  double sig_u = 1.0;
  double tail_scale = 1.0;
  arma::vec d_ar, s_ar, phi;
  std::vector<int> bnd_idx;
  std::vector<BoundaryRow> bnd_rows;
  arma::mat l_sub;

  void init(const arma::uvec& keep, double rho_in,
            const Rcpp::Nullable<Rcpp::NumericVector>& ar_phi,
            const Rcpp::Nullable<Rcpp::NumericVector>& ar_acf,
            double sig_u2,
            bool unit_marginal_in = false) {
    rho = rho_in;
    is_ar_p = false;
    use_dense_fallback = false;
    unit_marginal = unit_marginal_in;
    tail_scale = 1.0;
    bnd_idx.clear();
    bnd_rows.clear();
    if (ar_phi.isNotNull()) {
      phi = Rcpp::as<arma::vec>(ar_phi.get());
      if (phi.n_elem > 0) is_ar_p = true;
    }
    if (!is_ar_p) {
      hrf_whiten_coefs(keep, rho, d_ar, s_ar, unit_marginal);
      return;
    }
    if (ar_acf.isNull()) Rcpp::stop("ar_acf is required when ar_phi is supplied");
    arma::vec acf = Rcpp::as<arma::vec>(ar_acf.get());
    if (!phi.is_finite() || !acf.is_finite() || acf.n_elem < 1 || !std::isfinite(sig_u2) || sig_u2 <= 0)
      Rcpp::stop("Invalid AR(p) whitening parameters");
    sig_u = std::sqrt(std::max(1e-12, sig_u2));
    tail_scale = unit_marginal ? (1.0 / sig_u) : 1.0;
    const arma::uword m = keep.n_elem, p = phi.n_elem;
    if (m == 0) return;

    bnd_idx.assign(m, -1);
    bool markov_ok = true;
    for (arma::uword i = 0; i < m; ++i) {
      if (i >= p && (keep(i) - keep(i - p)) == p) {
        continue;
      }
      arma::uword b0 = 0;
      if (p > 0 && i >= p) {
        for (int j = static_cast<int>(i) - 1; j >= static_cast<int>(p) - 1; --j) {
          arma::uword uj = static_cast<arma::uword>(j);
          if ((keep(uj) - keep(uj - p + 1)) == (p - 1)) {
            b0 = uj - p + 1;
            break;
          }
        }
      }
      arma::uword s = i - b0 + 1;
      BoundaryRow br;
      br.b0 = b0;
      if (s == 1) {
        if (acf(0) <= 0.0) {
          markov_ok = false;
          break;
        }
        br.w_unit = arma::vec({1.0 / std::sqrt(acf(0))});
      } else {
        arma::mat g_loc(s, s);
        for (arma::uword a = 0; a < s; ++a) {
          g_loc(a, a) = acf(0);
          for (arma::uword b = a + 1; b < s; ++b) {
            arma::uword lag = keep(b0 + b) - keep(b0 + a);
            double v = (lag < acf.n_elem) ? acf(lag) : 0.0;
            g_loc(a, b) = v;
            g_loc(b, a) = v;
          }
        }
        arma::mat l_loc;
        if (!arma::chol(l_loc, g_loc, "lower")) {
          markov_ok = false;
          break;
        }
        arma::vec e_last(s, arma::fill::zeros);
        e_last(s - 1) = 1.0;
        if (!arma::solve(br.w_unit, arma::trimatu(l_loc.t()), e_last) ||
            !br.w_unit.is_finite() || br.w_unit(s - 1) <= 0.0) {
          markov_ok = false;
          break;
        }
      }
      bnd_idx[i] = static_cast<int>(bnd_rows.size());
      bnd_rows.push_back(std::move(br));
    }
    if (markov_ok) return;

    use_dense_fallback = true;
    arma::mat sigma_sub(m, m);
    for (arma::uword a = 0; a < m; ++a) {
      sigma_sub(a, a) = acf(0);
      for (arma::uword b = a + 1; b < m; ++b) {
        arma::uword lag = keep(b) - keep(a);
        double v = (lag < acf.n_elem) ? acf(lag) : 0.0;
        sigma_sub(a, b) = v;
        sigma_sub(b, a) = v;
      }
    }
    if (!arma::chol(l_sub, sigma_sub, "lower")) {
      sigma_sub.diag() += 1e-6;
      if (!arma::chol(l_sub, sigma_sub, "lower"))
        Rcpp::stop("AR(p) subset covariance Cholesky failed");
    }
  }

  void apply(arma::mat& x) const {
    if (x.n_cols == 0 || x.n_rows == 0) return;
    if (!is_ar_p) {
      if (x.n_rows == 1) {
        if (rho != 0) x.row(0) *= d_ar(0);
      } else {
        hrf_whiten_inplace(x, d_ar, s_ar, rho);
      }
      return;
    }
    const double head_scale = unit_marginal ? 1.0 : sig_u;
    if (use_dense_fallback) {
      x = head_scale * arma::solve(arma::trimatl(l_sub), x, arma::solve_opts::fast);
      return;
    }
    const arma::uword m = x.n_rows, p = phi.n_elem;
    const double* pp = phi.memptr();
    const double sc = tail_scale;
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      double* col = x.colptr(j);
      for (arma::uword ii = m; ii-- > 0; ) {
        int bi = bnd_idx[ii];
        if (bi < 0) {
          double acc = col[ii];
          for (arma::uword k = 0; k < p; ++k) {
            acc -= pp[k] * col[ii - 1 - k];
          }
          col[ii] = acc * sc;
        } else {
          const BoundaryRow& br = bnd_rows[static_cast<size_t>(bi)];
          const double* wp = br.w_unit.memptr();
          const arma::uword s = br.w_unit.n_elem;
          const arma::uword b0 = br.b0;
          double acc = 0.0;
          for (arma::uword a = 0; a < s; ++a) {
            acc += wp[a] * col[b0 + a];
          }
          col[ii] = acc * head_scale;
        }
      }
    }
  }

  // Compute U = L^T * bt = W_unit^{-T} * bt in-place so U^T * U = B * Sigma_sub * B^T
  // in O(m * p * M) time and O(1) extra memory (requires unit_marginal == true on AR(1)).
  void apply_inv_transpose(arma::mat& bt) const {
    if (bt.n_cols == 0 || bt.n_rows == 0) return;
    const arma::uword m = bt.n_rows;
    if (!is_ar_p) {
      if (rho == 0) return;
      const double* dp = d_ar.memptr();
      const double* sp = s_ar.memptr();
      for (arma::uword j = 0; j < bt.n_cols; ++j) {
        double* col = bt.colptr(j);
        for (arma::uword ii = m - 1; ii >= 1; --ii) {
          double u_i = col[ii] / dp[ii];
          col[ii] = u_i;
          col[ii - 1] -= sp[ii] * u_i;
        }
        col[0] /= dp[0];
      }
      return;
    }
    if (use_dense_fallback) {
      bt = l_sub.t() * bt;
      return;
    }
    const arma::uword p = phi.n_elem;
    const double* pp = phi.memptr();
    for (arma::uword j = 0; j < bt.n_cols; ++j) {
      double* col = bt.colptr(j);
      for (arma::uword ii = m; ii-- > 0; ) {
        int bi = bnd_idx[ii];
        if (bi < 0) {
          double v_i = col[ii];
          col[ii] = sig_u * v_i;
          for (arma::uword k = 0; k < p; ++k) {
            col[ii - 1 - k] += pp[k] * v_i;
          }
        } else {
          const BoundaryRow& br = bnd_rows[static_cast<size_t>(bi)];
          const double* wp = br.w_unit.memptr();
          const arma::uword s = br.w_unit.n_elem;
          const arma::uword b0 = br.b0;
          double u_i = col[ii] / wp[s - 1];
          col[ii] = u_i;
          for (arma::uword a = 0; a + 1 < s; ++a) {
            col[b0 + a] -= wp[a] * u_i;
          }
        }
      }
    }
  }
};

arma::uvec hrf_obs_rows(const Rcpp::IntegerVector& rows) {
  if (rows.size() < 1) return arma::uvec();
  arma::uvec out(rows.size());
  for (int i = 0; i < rows.size(); ++i) {
    if (rows[i] == NA_INTEGER || rows[i] < 1 || (i > 0 && rows[i] <= rows[i - 1]))
      Rcpp::stop("Observation rows must be positive and strictly increasing");
    out(i) = static_cast<arma::uword>(rows[i] - 1);
  }
  return out;
}

arma::mat hrf_extract_rows(const arma::mat& src_mat, const arma::uvec& keep,
                           arma::uword start, arma::uword end) {
  arma::uword m = keep.n_elem, b = end - start + 1;
  arma::mat y(m, b);
  if (m == src_mat.n_rows) {
    std::memcpy(y.memptr(), src_mat.colptr(start), m * b * sizeof(double));
  } else {
    const arma::uword* kp = keep.memptr();
    for (arma::uword j = 0; j < b; ++j) {
      const double* src = src_mat.colptr(start + j);
      double* dst = y.colptr(j);
      for (arma::uword i = 0; i < m; ++i) dst[i] = src[kp[i]];
    }
  }
  return y;
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
    bool evidence = false, int fallback_index = 1,
    Rcpp::Nullable<Rcpp::NumericVector> ar_phi = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> ar_acf = R_NilValue,
    double sig_u2 = 1.0) {
  if (bold.n_rows != nuisance.n_rows || bold.n_cols < 1 || rows.size() < 3 ||
      rows.size() > bold.n_rows || batch_size < 1 || designs.size() < 1 ||
      designs.size() > std::numeric_limits<int>::max() || rank_tol <= 0)
    Rcpp::stop("Invalid HRF scoring dimensions");
  if (!std::isfinite(rho) || std::abs(rho) >= 1 || fallback_index < 1 || fallback_index > designs.size())
    Rcpp::stop("Invalid AR coefficient or fallback index");
  arma::uvec keep = hrf_rows(rows, bold.n_rows);
  ArWhitener whitener;
  whitener.init(keep, rho, ar_phi, ar_acf, sig_u2);
  arma::mat z = nuisance.rows(keep);
  if (!z.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
  whitener.apply(z);
  arma::mat qn = hrf_basis(z, rank_tol);
  int base_df = keep.n_elem - qn.n_cols;
  const int n_cand = static_cast<int>(designs.size());
  std::vector<arma::mat> qt(n_cand);
  Rcpp::IntegerVector ranks(n_cand);
  for (int k = 0; k < n_cand; ++k) {
    arma::mat x = Rcpp::as<arma::mat>(designs[k]);
    if (x.n_rows != bold.n_rows || !x.is_finite()) Rcpp::stop("Invalid event design");
    x = x.rows(keep);
    whitener.apply(x);
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
  Rcpp::NumericMatrix score(n_cand, bold.n_cols);
  std::fill(score.begin(), score.end(), NA_REAL);
  int evidence_rows = evidence ? n_cand : 0, evidence_cols = evidence ? bold.n_cols : 0;
  Rcpp::NumericMatrix task_p(evidence_rows, evidence_cols), fallback_p(evidence_rows, evidence_cols),
    incompatible_p(evidence_rows, evidence_cols);
  std::fill(task_p.begin(), task_p.end(), NA_REAL);
  std::fill(fallback_p.begin(), fallback_p.end(), NA_REAL);
  std::fill(incompatible_p.begin(), incompatible_p.end(), NA_REAL);

  const bool exact_all_pairs = (n_cand <= 18);
  int fb = fallback_index - 1;
  struct Pair { int i, j, extra_rank, df; arma::mat q; };
  std::vector<Pair> pairs;
  std::vector<int> pair_lookup;
  std::vector<int> pair_seen;
  std::vector<int> global_anchors;

  auto get_pair = [&](int a, int b) -> int {
    if (a == b) return -2;
    int i = std::min(a, b), j = std::max(a, b);
    int key = i * n_cand + j;
    int st = pair_lookup[key];
    if (st != -1) return st;
    arma::mat extra = qt[j] - qt[i] * (qt[i].t() * qt[j]);
    for (arma::uword c = 0; c < extra.n_cols; ++c)
      if (arma::norm(extra.col(c)) <= rank_tol) extra.col(c).zeros();
    arma::mat q = hrf_basis(extra, rank_tol);
    int d = q.n_cols;
    if (d && base_df - ranks[i] - d >= 2) {
      int idx = static_cast<int>(pairs.size());
      pairs.push_back(Pair{i, j, d, base_df - ranks[i] - d, std::move(q)});
      pair_seen.push_back(-1);
      pair_lookup[key] = idx;
      return idx;
    }
    pair_lookup[key] = -2;
    return -2;
  };

  if (comparable && evidence) {
    if (exact_all_pairs) {
      for (int i = 0; i < n_cand; ++i)
        for (int j = i + 1; j < n_cand; ++j) {
          // Normalize only genuinely nonzero residual columns, so equal/subspace-
          // equivalent candidates never acquire rank from roundoff.
          arma::mat extra = qt[j] - qt[i] * (qt[i].t() * qt[j]);
          for (arma::uword c = 0; c < extra.n_cols; ++c)
            if (arma::norm(extra.col(c)) <= rank_tol) extra.col(c).zeros();
          arma::mat q = hrf_basis(extra, rank_tol);
          int d = q.n_cols;
          if (d && base_df - ranks[i] - d >= 2)
            pairs.push_back(Pair{i, j, d, base_df - ranks[i] - d, std::move(q)});
        }
    } else {
      pair_lookup.assign(static_cast<size_t>(n_cand) * static_cast<size_t>(n_cand), -1);
      global_anchors.push_back(fb);
      arma::vec max_ov(n_cand, arma::fill::zeros);
      max_ov(fb) = 1.0;
      double inv_rk = 1.0 / std::max(1, static_cast<int>(ranks[0]));
      for (int k = 0; k < n_cand; ++k) {
        if (k == fb) continue;
        max_ov(k) = arma::accu(arma::square(qt[fb].t() * qt[k])) * inv_rk;
      }
      int n_medoids = std::min(8, n_cand - 1);
      for (int m_idx = 0; m_idx < n_medoids; ++m_idx) {
        arma::uword next_k = max_ov.index_min();
        if (max_ov(next_k) >= 1.0 - 1e-10) break;
        int ak = static_cast<int>(next_k);
        global_anchors.push_back(ak);
        max_ov(ak) = 1.0;
        for (int k = 0; k < n_cand; ++k) {
          if (max_ov(k) >= 1.0 - 1e-10) continue;
          double ov = arma::accu(arma::square(qt[ak].t() * qt[k])) * inv_rk;
          if (ov > max_ov(k)) max_ov(k) = ov;
        }
      }
    }
  }
  struct FallbackAngle { int j, df; arma::mat pin_w, u_orth; };
  std::vector<FallbackAngle> fb_angles;
  if (comparable && evidence && ranks[0] > 3 && fb >= 0 && fb < n_cand &&
      base_df - ranks[fb] - 1 >= 2) {
    const arma::mat& Qi = qt[fb];
    int df_fb = base_df - ranks[fb] - 1;
    for (int j = 0; j < n_cand; ++j) {
      if (j == fb) continue;
      const arma::mat& Qj = qt[j];
      arma::mat U_c, V_c;
      arma::vec s_c;
      if (!arma::svd(U_c, s_c, V_c, Qi.t() * Qj)) continue;
      arma::vec cos_t = arma::clamp(s_c, 0.0, 1.0);
      arma::vec sin_t = arma::sqrt(arma::clamp(1.0 - arma::square(cos_t), 0.0, 1.0));
      double max_sin = sin_t.max();
      if (max_sin <= 1e-6) continue;
      arma::uvec keep_idx = arma::find(sin_t >= 0.5 * max_sin && cos_t > 0.05);
      if (keep_idx.is_empty()) {
        keep_idx = arma::uvec{ sin_t.index_max() };
      }
      arma::mat Uc = U_c.cols(keep_idx);
      arma::mat Vc = V_c.cols(keep_idx);
      arma::vec cos_k = cos_t.elem(keep_idx);
      arma::vec sin_k = sin_t.elem(keep_idx);
      arma::vec tan_k = sin_k / arma::clamp(cos_k, 0.05, 1.0);
      arma::mat Pin = Qi * Uc;
      arma::mat Uorth = Qj * Vc - Pin * arma::diagmat(cos_k);
      Uorth.each_row() /= sin_k.t();
      arma::mat Pin_w = Pin * arma::diagmat(tan_k);
      fb_angles.push_back(FallbackAngle{j, df_fb, std::move(Pin_w), std::move(Uorth)});
    }
  }
  Rcpp::CharacterVector status(bold.n_cols, "unidentifiable_design");
  int batch_id = 0;
  if (comparable) for (arma::uword start = 0; start < bold.n_cols; start += batch_size, ++batch_id) {
    Rcpp::checkUserInterrupt();
    arma::uword end = std::min(bold.n_cols - 1, start + batch_size - 1);
    arma::mat y = hrf_extract_rows(bold, keep, start, end);
    std::vector<uint8_t> ok_vox(y.n_cols, 0);
    for (arma::uword j = 0; j < y.n_cols; ++j) {
      if (!y.col(j).is_finite()) {
        status[start + j] = "nonfinite";
        y.col(j).zeros();
      } else {
        y.col(j) -= arma::mean(y.col(j));
        double norm = arma::norm(y.col(j));
        if (norm > 0) {
          status[start + j] = "ok";
          ok_vox[j] = 1;
          y.col(j) /= norm;
        } else {
          status[start + j] = "constant";
        }
      }
    }
    whitener.apply(y);
    if (qn.n_cols) y -= qn * (qn.t() * y);
    arma::rowvec sst = arma::sum(arma::square(y), 0);
    for (arma::uword j = 0; j < y.n_cols; ++j) {
      if (ok_vox[j] && sst(j) <= rank_tol * rank_tol) {
        status[start + j] = "nuisance_only";
        ok_vox[j] = 0;
      }
      if (sst(j) > rank_tol * rank_tol) y.col(j) /= std::sqrt(sst(j));
    }
    arma::mat all_explained(n_cand, y.n_cols, arma::fill::zeros);
    for (int k = 0; k < n_cand; ++k) {
      arma::rowvec explained = arma::sum(arma::square(qt[k].t() * y), 0);
      for (arma::uword j = 0; j < y.n_cols; ++j) if (ok_vox[j]) {
        double r2 = std::max(0.0, std::min(1.0, explained(j)));
        score(k, start + j) = 1 - (1 - r2) * base_df / (base_df - ranks[k]);
        all_explained(k, j) = r2;
        if (evidence) {
          task_p(k, start + j) = hrf_f_tail(r2, ranks[k], 1 - r2, base_df - ranks[k]);
          fallback_p(k, start + j) = incompatible_p(k, start + j) = 1;
        }
      }
    }
    auto eval_pair = [&](const Pair& pair) {
      arma::rowvec added = arma::sum(arma::square(pair.q.t() * y), 0);
      for (arma::uword j = 0; j < y.n_cols; ++j)
        if (ok_vox[j]) {
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
          if (pair.i == fb) fallback_p(pair.j, start + j) = pi;
          if (pair.j == fb) fallback_p(pair.i, start + j) = pj;
        }
    };
    if (evidence) {
      if (exact_all_pairs) {
        for (const Pair& pair : pairs) eval_pair(pair);
      } else {
        std::vector<int> win_cnt(n_cand, 0);
        std::vector<double> win_max_r2(n_cand, -1.0);
        for (arma::uword j = 0; j < y.n_cols; ++j) if (ok_vox[j]) {
          int best_k = static_cast<int>(all_explained.col(j).index_max());
          win_cnt[best_k] += 1;
          if (all_explained(best_k, j) > win_max_r2[best_k])
            win_max_r2[best_k] = all_explained(best_k, j);
        }
        std::vector<int> winners;
        for (int k = 0; k < n_cand; ++k) if (win_cnt[k] > 0) winners.push_back(k);
        if (winners.size() > 10) {
          std::sort(winners.begin(), winners.end(), [&](int a, int b) {
            if (win_max_r2[a] != win_max_r2[b]) return win_max_r2[a] > win_max_r2[b];
            if (win_cnt[a] != win_cnt[b]) return win_cnt[a] > win_cnt[b];
            return a < b;
          });
          winners.resize(10);
        }
        std::vector<uint8_t> is_anchor(n_cand, 0);
        std::vector<int> batch_anchors;
        batch_anchors.reserve(global_anchors.size() + winners.size());
        for (int a : global_anchors) if (!is_anchor[a]) { is_anchor[a] = 1; batch_anchors.push_back(a); }
        for (int a : winners) if (!is_anchor[a]) { is_anchor[a] = 1; batch_anchors.push_back(a); }

        std::vector<int> active_pids;
        active_pids.reserve(batch_anchors.size() * static_cast<size_t>(n_cand));
        for (int a : batch_anchors) {
          for (int k = 0; k < n_cand; ++k) {
            if (k == a) continue;
            int pid = get_pair(a, k);
            if (pid >= 0 && pair_seen[pid] != batch_id) {
              pair_seen[pid] = batch_id;
              active_pids.push_back(pid);
            }
          }
        }
        for (int pid : active_pids) eval_pair(pairs[pid]);
      }
    }
    if (evidence && !fb_angles.empty()) {
      for (const FallbackAngle& fa : fb_angles) {
        arma::mat w = fa.pin_w.t() * y;
        arma::mat b = fa.u_orth.t() * y;
        for (arma::uword j = 0; j < y.n_cols; ++j) {
          if (ok_vox[j]) {
            double w_norm = arma::norm(w.col(j));
            if (w_norm <= 1e-12) {
              fallback_p(fa.j, start + j) = 1.0;
              continue;
            }
            double u_proj = arma::dot(w.col(j), b.col(j)) / w_norm;
            double extra_ss = u_proj * u_proj;
            double ei = all_explained(fb, j);
            double both = std::min(1.0, std::max(ei, ei + extra_ss));
            fallback_p(fa.j, start + j) = hrf_f_tail(both - ei, 1, 1.0 - both, fa.df);
          }
        }
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
  struct RhoFactor {
    double rho;
    arma::vec d_ar, s_ar;
    arma::mat q;
    double log_term;
  };
  std::vector<RhoFactor> factors(rho_grid.n_elem);
  for (arma::uword k = 0; k < rho_grid.n_elem; ++k) {
    Rcpp::checkUserInterrupt();
    double rho = rho_grid(k);
    if (!std::isfinite(rho) || std::abs(rho) >= 1) Rcpp::stop("Invalid noise grid");
    RhoFactor& f = factors[k];
    f.rho = rho;
    hrf_whiten_coefs(keep, rho, f.d_ar, f.s_ar);
    arma::mat wx = basis, v;
    hrf_whiten_inplace(wx, f.d_ar, f.s_ar, rho);
    arma::vec s;
    if (!arma::svd_econ(f.q, s, v, wx)) Rcpp::stop("Noise design SVD failed");
    double logdet_design = 2 * arma::sum(arma::log(s));
    double logdet_w = .5 * std::log1p(-rho * rho);
    for (arma::uword i = 1; i < keep.n_elem; ++i) {
      double a = std::pow(rho, (double)(keep(i) - keep(i - 1)));
      logdet_w += .5 * (std::log1p(-rho * rho) - std::log1p(-a * a));
    }
    f.log_term = logdet_design - 2 * logdet_w;
  }
  x.reset();
  basis.reset();
  arma::mat y;
  for (arma::uword start = 0; start < bold.n_cols; start += batch_size) {
    Rcpp::checkUserInterrupt();
    arma::uword end = std::min(bold.n_cols - 1, start + batch_size - 1);
    arma::mat y_base = hrf_extract_rows(bold, keep, start, end);
    arma::uvec valid(y_base.n_cols, arma::fill::ones);
    for (arma::uword j = 0; j < y_base.n_cols; ++j) {
      if (!y_base.col(j).is_finite()) { valid(j) = 0; y_base.col(j).zeros(); continue; }
      y_base.col(j) -= arma::mean(y_base.col(j));
      double norm = arma::norm(y_base.col(j));
      if (norm <= 0) valid(j) = 0;
      else y_base.col(j) /= norm;
    }
    for (arma::uword k = 0; k < rho_grid.n_elem; ++k) {
      const RhoFactor& f = factors[k];
      y = y_base;
      hrf_whiten_inplace(y, f.d_ar, f.s_ar, f.rho);
      arma::rowvec total = arma::sum(arma::square(y), 0), rss = total;
      if (f.q.n_cols) {
        arma::mat coefficients = f.q.t() * y;
        rss -= arma::sum(arma::square(coefficients), 0);
        // Avoid a second large BLAS multiplication for ordinary-SNR voxels.
        // Explicit residuals protect near-perfect fits from cancellation in
        // ||y||^2 - ||Q'y||^2 (including negative roundoff differences).
        for (arma::uword j = 0; j < y.n_cols; ++j)
          if (valid(j) && rss(j) < 1e-4 * total(j)) {
            arma::vec residual = y.col(j) - f.q * coefficients.col(j);
            rss(j) = arma::dot(residual, residual);
          }
      }
      for (arma::uword j = 0; j < y.n_cols; ++j) if (valid(j) && rss(j) > rank_tol * rank_tol)
        profile(k, start + j) = df * std::log(rss(j) / df) + f.log_term;
    }
  }
  return profile;
}

// Exact masked GLS inverse, jointly fitting unpenalized nuisance coefficients,
// with optional Diedrichsen-style / Woodbury rWLS observation weights `obs_weights`
// in (0, 1] on `rows` (`weight_space = "woodbury"` models additive transient variance
// Sigma_sub + diag((1 - w) / w); `weight_space = "innovation"` scales whitened rows
// by sqrt(w)). Return only its action on measured rows; optional padded zeros are fixed.
// [[Rcpp::export]]
arma::mat hrf_inverse_operator_cpp(int n_time, const arma::vec& kernel,
    const Rcpp::IntegerVector& rows, double rho, const arma::mat& nuisance,
    double lambda, const std::string& penalty, double ridge_floor, int tail_padding = 0,
    Rcpp::Nullable<Rcpp::NumericVector> ar_phi = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> ar_acf = R_NilValue,
    double sig_u2 = 1.0,
    Rcpp::Nullable<Rcpp::NumericVector> obs_weights = R_NilValue,
    const std::string& weight_space = "woodbury") {
  arma::uvec keep = hrf_rows(rows, n_time);
  if (kernel.n_elem < 2 || !kernel.is_finite() || nuisance.n_rows != (arma::uword)n_time ||
      !std::isfinite(rho) || std::abs(rho) >= 1 || lambda < 0 || ridge_floor < 0 || tail_padding < 0)
    Rcpp::stop("Invalid masked inverse parameters");
  if (weight_space != "woodbury" && weight_space != "innovation")
    Rcpp::stop("Unknown weight_space in hrf_inverse_operator_cpp");
  double pen_w0 = 0.0, pen_w1 = 0.0, pen_w2 = 0.0;
  arma::uword pen_band = 0;
  if (penalty == "ridge") {
    pen_w0 = 1.0; pen_band = 0;
  } else if (penalty == "diff1") {
    pen_w1 = 1.0; pen_band = 1;
  } else if (penalty == "diff2") {
    pen_w2 = 1.0; pen_band = 2;
  } else if (penalty == "matern1") {
    pen_w0 = 0.10; pen_w1 = 0.90; pen_band = 1;
  } else if (penalty == "sobolev") {
    pen_w0 = 0.05; pen_w1 = 0.45; pen_w2 = 0.50; pen_band = 2;
  } else {
    Rcpp::stop("Unknown inverse penalty");
  }

  arma::uword measured = keep.n_elem, k_len = kernel.n_elem, history = k_len - 1;
  if (tail_padding) keep = arma::join_cols(keep,
    arma::regspace<arma::uvec>(n_time, n_time + tail_padding - 1));
  arma::uword m = keep.n_elem, latent = n_time + tail_padding + history;

  arma::vec w_obs;
  bool has_obs_w = false;
  if (obs_weights.isNotNull()) {
    w_obs = Rcpp::as<arma::vec>(obs_weights.get());
    if (w_obs.n_elem != measured || !w_obs.is_finite() || arma::any(w_obs <= 0.0))
      Rcpp::stop("obs_weights must be positive, finite, and match length(rows)");
    for (arma::uword i = 0; i < measured; ++i) {
      if (w_obs(i) > 1.0) w_obs(i) = 1.0;
      if (std::abs(w_obs(i) - 1.0) > 1e-8) has_obs_w = true;
    }
  }

  arma::vec k_rev = arma::reverse(kernel);
  const double* kp = k_rev.memptr();

  // Cascaded AR(p) OR rWLS weighted path: fast in-place AR(p)/AR(1) whitener +
  // optional low-rank Woodbury square-root update + BLAS DGEMM + LAPACK Cholesky solve.
  ArWhitener whitener;
  const bool use_woodbury = has_obs_w && (weight_space == "woodbury");
  whitener.init(keep, rho, ar_phi, ar_acf, sig_u2, use_woodbury);
  if (whitener.is_ar_p || has_obs_w) {
    arma::mat u_a;
    arma::vec sh_a;
    arma::vec sqrt_w_all;
    double s0_wood = 1.0;
    if (use_woodbury) {
      s0_wood = whitener.is_ar_p
        ? std::sqrt(std::max(sig_u2, 1e-8))
        : (rho == 0.0 ? 1.0 : std::sqrt(std::max(1.0 - rho * rho, 1e-8)));
      std::vector<arma::uword> sp_idx;
      sp_idx.reserve(measured);
      for (arma::uword i = 0; i < measured; ++i) {
        if (w_obs(i) < 1.0 - 1e-8) sp_idx.push_back(i);
      }
      if (!sp_idx.empty()) {
        arma::mat a_mat(m, sp_idx.size(), arma::fill::zeros);
        for (size_t c = 0; c < sp_idx.size(); ++c) {
          arma::uword ii = sp_idx[c];
          double w_ii = std::max(w_obs(ii), 1e-8);
          a_mat(ii, c) = std::sqrt((1.0 - w_ii) / w_ii);
        }
        whitener.apply(a_mat);
        arma::vec d_a;
        arma::mat v_a;
        bool ok_svd = arma::svd_econ(u_a, d_a, v_a, a_mat, "left");
        if (!ok_svd) Rcpp::stop("SVD failed in Woodbury rWLS whitener");
        sh_a = 1.0 - 1.0 / arma::sqrt(1.0 + arma::square(d_a));
      }
    } else if (has_obs_w) {
      sqrt_w_all.ones(m);
      sqrt_w_all.head(measured) = arma::sqrt(arma::clamp(w_obs, 1e-8, 1.0));
    }

    auto apply_eff_whiten = [&](arma::mat& mat_x) {
      whitener.apply(mat_x);
      if (use_woodbury) {
        if (u_a.n_cols > 0) {
          arma::mat utx = u_a.t() * mat_x;
          utx.each_col() %= sh_a;
          mat_x -= u_a * utx;
        }
        if (s0_wood != 1.0) mat_x *= s0_wood;
      } else if (has_obs_w) {
        mat_x.each_col() %= sqrt_w_all;
      }
    };

    arma::mat q;
    if (nuisance.n_cols > 0) {
      arma::mat z(m, nuisance.n_cols, arma::fill::zeros);
      z.rows(0, measured - 1) = nuisance.rows(keep.head(measured));
      if (!z.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
      apply_eff_whiten(z);
      q = hrf_basis(z, 1e-8);
    } else {
      arma::mat z_check = nuisance.rows(keep.head(measured));
      if (!z_check.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
    }
    arma::mat h(m, latent, arma::fill::zeros);
    for (arma::uword i = 0; i < m; ++i) {
      arma::uword r0 = keep(i);
      for (arma::uword u = 0; u < k_len; ++u) h(i, r0 + u) = kp[u];
    }
    apply_eff_whiten(h);
    arma::mat w_meas(m, measured, arma::fill::zeros);
    w_meas.rows(0, measured - 1).eye();
    apply_eff_whiten(w_meas);
    if (q.n_cols > 0) {
      h -= q * (q.t() * h);
      w_meas -= q * (q.t() * w_meas);
    }
    arma::mat lhs = h.t() * h;
    arma::mat rhs = h.t() * w_meas;
    if (pen_w0 > 0.0) lhs.diag() += lambda * pen_w0;
    if (pen_w1 > 0.0) {
      const double d1[] = {-1, 1};
      for (arma::uword i = 0; i + 1 < latent; ++i)
        for (int a = 0; a <= 1; ++a)
          for (int b = a; b <= 1; ++b)
            lhs(i + b, i + a) += (lambda * pen_w1) * d1[a] * d1[b];
    }
    if (pen_w2 > 0.0) {
      const double d2[] = {1, -2, 1};
      for (arma::uword i = 0; i + 2 < latent; ++i)
        for (int a = 0; a <= 2; ++a)
          for (int b = a; b <= 2; ++b)
            lhs(i + b, i + a) += (lambda * pen_w2) * d2[a] * d2[b];
    }
    lhs.diag() += ridge_floor;

    arma::blas_int n_lat = static_cast<arma::blas_int>(latent);
    arma::blas_int n_rhs = static_cast<arma::blas_int>(measured);
    arma::blas_int info = 0;
    char uplo = 'L';
    arma::lapack::potrf(&uplo, &n_lat, lhs.memptr(), &n_lat, &info);
    if (info != 0) Rcpp::stop("Masked GLS inverse failed; increase ridge_floor");
    arma::lapack::potrs(&uplo, &n_lat, &n_rhs, lhs.memptr(), &n_lat, rhs.memptr(), &n_lat, &info);
    if (info != 0) Rcpp::stop("Masked GLS inverse failed; increase ridge_floor");
    return rhs.rows(history, history + n_time - 1);
  }

  // Bidiagonal AR(1) whitening coefficients on the observation clock (no m x m matrix)
  const arma::vec& d = whitener.d_ar;
  const arma::vec& s = whitener.s_ar;

  arma::mat q;
  if (nuisance.n_cols > 0) {
    arma::mat z(m, nuisance.n_cols, arma::fill::zeros);
    z.rows(0, measured - 1) = nuisance.rows(keep.head(measured));
    if (!z.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
    hrf_whiten_inplace(z, d, s, rho);
    q = hrf_basis(z, 1e-8);
  } else {
    arma::mat z_check = nuisance.rows(keep.head(measured));
    if (!z_check.is_finite()) Rcpp::stop("Nuisance values must be finite on usable rows");
  }

  // Symmetric tridiagonal entries of M = W^T W
  arma::vec m_diag(m), m_sub(m, arma::fill::zeros);
  arma::uword max_gap = 0;
  for (arma::uword i = 0; i < m; ++i) {
    m_diag(i) = d(i) * d(i) + (i + 1 < m ? s(i + 1) * s(i + 1) : 0.0);
    if (i > 0) {
      m_sub(i) = d(i) * s(i);
      if (rho != 0) max_gap = std::max(max_gap, keep(i) - keep(i - 1));
    }
  }

  arma::mat g = k_rev * k_rev.t();

  // Fast banded Cholesky + L1-tiled SIMD solve when no dense nuisance projection is present
  arma::uword kd = std::max(max_gap + history, pen_band);
  if (q.n_cols == 0 && kd < latent / 3) {
    arma::mat ab(kd + 1, latent, arma::fill::zeros);
    for (arma::uword i = 0; i < m; ++i) {
      arma::uword r0 = keep(i);
      double md = m_diag(i);
      for (arma::uword u = 0; u < k_len; ++u) {
        arma::uword col_idx = r0 + u;
        for (arma::uword v = u; v < k_len; ++v) {
          ab(v - u, col_idx) += md * g(v, u);
        }
      }
      if (i > 0 && m_sub(i) != 0) {
        arma::uword rp = keep(i - 1);
        double ms = m_sub(i);
        for (arma::uword u = 0; u < k_len; ++u) {
          arma::uword r_curr = r0 + u;
          for (arma::uword v = 0; v < k_len; ++v) {
            arma::uword r_prev = rp + v;
            double val = ms * g(u, v);
            if (r_curr >= r_prev) ab(r_curr - r_prev, r_prev) += val;
            if (r_prev >= r_curr) ab(r_prev - r_curr, r_curr) += val;
          }
        }
      }
    }
    if (pen_w0 > 0.0) ab.row(0) += lambda * pen_w0;
    if (pen_w1 > 0.0) {
      const double d1[] = {-1, 1};
      for (arma::uword i = 0; i + 1 < latent; ++i)
        for (int a = 0; a <= 1; ++a)
          for (int b = a; b <= 1; ++b)
            ab(b - a, i + a) += (lambda * pen_w1) * d1[a] * d1[b];
    }
    if (pen_w2 > 0.0) {
      const double d2[] = {1, -2, 1};
      for (arma::uword i = 0; i + 2 < latent; ++i)
        for (int a = 0; a <= 2; ++a)
          for (int b = a; b <= 2; ++b)
            ab(b - a, i + a) += (lambda * pen_w2) * d2[a] * d2[b];
    }
    ab.row(0) += ridge_floor;

    arma::blas_int n_lat = static_cast<arma::blas_int>(latent);
    arma::blas_int kd_int = static_cast<arma::blas_int>(kd);
    arma::blas_int ldab = kd_int + 1;
    arma::blas_int info = 0;
    char uplo = 'L';
    arma::lapack::pbtrf(&uplo, &n_lat, &kd_int, ab.memptr(), &ldab, &info);
    if (info == 0) {
      for (arma::uword j = 0; j < latent; ++j) {
        ab(0, j) = 1.0 / ab(0, j);
      }
      arma::mat out(n_time, measured);
      constexpr arma::uword TILE = 32;
      arma::mat tile(TILE, latent);

      for (arma::uword c0 = 0; c0 < measured; c0 += TILE) {
        arma::uword c1 = std::min(measured, c0 + TILE);
        arma::uword b_cols = c1 - c0;
        arma::uword j_start = (c0 > 0 && m_sub(c0) != 0) ? keep(c0 - 1) : keep(c0);
        arma::uword clear_start = std::min(history, j_start);
        std::memset(tile.colptr(clear_start), 0, (latent - clear_start) * TILE * sizeof(double));

        for (arma::uword idx = 0; idx < b_cols; ++idx) {
          arma::uword j = c0 + idx;
          arma::uword r0 = keep(j);
          double md = m_diag(j);
          for (arma::uword u = 0; u < k_len; ++u) tile(idx, r0 + u) += md * kp[u];
          if (j > 0 && m_sub(j) != 0) {
            arma::uword rp = keep(j - 1);
            double ms = m_sub(j);
            for (arma::uword u = 0; u < k_len; ++u) tile(idx, rp + u) += ms * kp[u];
          }
          if (j + 1 < m && m_sub(j + 1) != 0) {
            arma::uword rn = keep(j + 1);
            double ms = m_sub(j + 1);
            for (arma::uword u = 0; u < k_len; ++u) tile(idx, rn + u) += ms * kp[u];
          }
        }

        for (arma::uword j = j_start; j < latent; ++j) {
          const double* l_col = ab.colptr(j);
          double inv_d = l_col[0];
          double* yj = tile.colptr(j);
          for (arma::uword i = 0; i < b_cols; ++i) yj[i] *= inv_d;
          arma::uword kmax = std::min(kd, latent - 1 - j);
          for (arma::uword k = 1; k <= kmax; ++k) {
            double lk = l_col[k];
            double* yjk = tile.colptr(j + k);
            for (arma::uword i = 0; i < b_cols; ++i) {
              yjk[i] -= lk * yj[i];
            }
          }
        }

        for (arma::uword jj = latent; jj > history; --jj) {
          arma::uword j = jj - 1;
          const double* l_col = ab.colptr(j);
          double inv_d = l_col[0];
          double* yj = tile.colptr(j);
          arma::uword kmax = std::min(kd, latent - 1 - j);
          for (arma::uword k = 1; k <= kmax; ++k) {
            double lk = l_col[k];
            const double* yjk = tile.colptr(j + k);
            for (arma::uword i = 0; i < b_cols; ++i) {
              yj[i] -= lk * yjk[i];
            }
          }
          for (arma::uword i = 0; i < b_cols; ++i) yj[i] *= inv_d;
        }

        for (arma::uword idx = 0; idx < b_cols; ++idx) {
          double* out_col = out.colptr(c0 + idx);
          for (arma::uword t = 0; t < (arma::uword)n_time; ++t) {
            out_col[t] = tile(idx, history + t);
          }
        }
      }
      return out;
    }
  }

  // General path (joint nuisance projection or wide-band gap) using stencil accumulation
  arma::mat lhs(latent, latent, arma::fill::zeros);
  for (arma::uword i = 0; i < m; ++i) {
    arma::uword r0 = keep(i);
    double md = m_diag(i);
    for (arma::uword u = 0; u < k_len; ++u) {
      double* col_ptr = lhs.colptr(r0 + u);
      for (arma::uword v = u; v < k_len; ++v) {
        col_ptr[r0 + v] += md * g(v, u);
      }
    }
    if (i > 0 && m_sub(i) != 0) {
      arma::uword rp = keep(i - 1);
      double ms = m_sub(i);
      for (arma::uword u = 0; u < k_len; ++u) {
        arma::uword r_curr = r0 + u;
        for (arma::uword v = 0; v < k_len; ++v) {
          arma::uword r_prev = rp + v;
          double val = ms * g(u, v);
          if (r_curr >= r_prev) lhs(r_curr, r_prev) += val;
          if (r_prev >= r_curr) lhs(r_prev, r_curr) += val;
        }
      }
    }
  }

  arma::mat rhs(latent, measured, arma::fill::zeros);
  for (arma::uword j = 0; j < measured; ++j) {
    double* col = rhs.colptr(j);
    arma::uword r0 = keep(j);
    double md = m_diag(j);
    for (arma::uword u = 0; u < k_len; ++u) col[r0 + u] += md * kp[u];
    if (j > 0 && m_sub(j) != 0) {
      arma::uword rp = keep(j - 1);
      double ms = m_sub(j);
      for (arma::uword u = 0; u < k_len; ++u) col[rp + u] += ms * kp[u];
    }
    if (j + 1 < m && m_sub(j + 1) != 0) {
      arma::uword rn = keep(j + 1);
      double ms = m_sub(j + 1);
      for (arma::uword u = 0; u < k_len; ++u) col[rn + u] += ms * kp[u];
    }
  }

  if (q.n_cols > 0) {
    arma::mat q_tilde(m, q.n_cols);
    for (arma::uword c = 0; c < q.n_cols; ++c) {
      const double* qc = q.colptr(c);
      double* qtc = q_tilde.colptr(c);
      for (arma::uword i = 0; i + 1 < m; ++i) qtc[i] = d(i) * qc[i] + s(i + 1) * qc[i + 1];
      qtc[m - 1] = d(m - 1) * qc[m - 1];
    }
    arma::vec uc(latent);
    for (arma::uword c = 0; c < q.n_cols; ++c) {
      uc.zeros();
      double* ucp = uc.memptr();
      const double* qtc = q_tilde.colptr(c);
      for (arma::uword i = 0; i < m; ++i) {
        double w_ic = qtc[i];
        arma::uword r0 = keep(i);
        for (arma::uword u = 0; u < k_len; ++u) ucp[r0 + u] += w_ic * kp[u];
      }
      for (arma::uword j = 0; j < latent; ++j) {
        double u_jc = ucp[j];
        double* lhs_col = lhs.colptr(j);
        for (arma::uword r = j; r < latent; ++r) {
          lhs_col[r] -= u_jc * ucp[r];
        }
      }
      for (arma::uword j = 0; j < measured; ++j) {
        double q_jc = qtc[j];
        double* rhs_col = rhs.colptr(j);
        for (arma::uword r = 0; r < latent; ++r) {
          rhs_col[r] -= q_jc * ucp[r];
        }
      }
    }
  }

  if (pen_w0 > 0.0) lhs.diag() += lambda * pen_w0;
  if (pen_w1 > 0.0) {
    const double d1[] = {-1, 1};
    for (arma::uword i = 0; i + 1 < latent; ++i)
      for (int a = 0; a <= 1; ++a)
        for (int b = a; b <= 1; ++b)
          lhs(i + b, i + a) += (lambda * pen_w1) * d1[a] * d1[b];
  }
  if (pen_w2 > 0.0) {
    const double d2[] = {1, -2, 1};
    for (arma::uword i = 0; i + 2 < latent; ++i)
      for (int a = 0; a <= 2; ++a)
        for (int b = a; b <= 2; ++b)
          lhs(i + b, i + a) += (lambda * pen_w2) * d2[a] * d2[b];
  }
  lhs.diag() += ridge_floor;

  arma::blas_int n_lat = static_cast<arma::blas_int>(latent);
  arma::blas_int n_rhs = static_cast<arma::blas_int>(measured);
  arma::blas_int info = 0;
  char uplo = 'L';
  arma::lapack::potrf(&uplo, &n_lat, lhs.memptr(), &n_lat, &info);
  if (info != 0) Rcpp::stop("Masked GLS inverse failed; increase ridge_floor");
  arma::lapack::potrs(&uplo, &n_lat, &n_rhs, lhs.memptr(), &n_lat, rhs.memptr(), &n_lat, &info);
  if (info != 0) Rcpp::stop("Masked GLS inverse failed; increase ridge_floor");
  return rhs.rows(history, history + n_time - 1);
}

// Apply exact AR(1) or cascaded AR(p) Cholesky whitening to a matrix (or list of
// matrices) already restricted to observed acquisition indices `rows` (m = length(rows)),
// without allocating a dense m x m whitening matrix in R.
// [[Rcpp::export]]
SEXP hrf_whiten_apply_cpp(SEXP target, const Rcpp::IntegerVector& rows,
    double rho = 0.0,
    Rcpp::Nullable<Rcpp::NumericVector> ar_phi = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> ar_acf = R_NilValue,
    double sig_u2 = 1.0,
    bool unit_marginal = false) {
  if (!std::isfinite(rho) || std::abs(rho) >= 1.0)
    Rcpp::stop("Invalid AR(1) coefficient in hrf_whiten_apply_cpp");
  arma::uvec keep = hrf_obs_rows(rows);
  const arma::uword m = keep.n_elem;
  ArWhitener whitener;
  if (m >= 1) {
    whitener.init(keep, rho, ar_phi, ar_acf, sig_u2, unit_marginal);
  }
  if (Rf_isNewList(target)) {
    Rcpp::List in_list(target);
    int n_el = in_list.size();
    Rcpp::List out_list(n_el);
    out_list.attr("names") = in_list.attr("names");
    for (int idx = 0; idx < n_el; ++idx) {
      if (Rf_isNull(in_list[idx])) {
        out_list[idx] = R_NilValue;
        continue;
      }
      arma::mat mat_i = Rcpp::as<arma::mat>(in_list[idx]);
      if (mat_i.n_rows != m)
        Rcpp::stop("Matrix row count in hrf_whiten_apply_cpp list must match length(rows)");
      if (m >= 1 && mat_i.n_cols > 0) whitener.apply(mat_i);
      out_list[idx] = mat_i;
    }
    return out_list;
  }
  arma::mat mat_x = Rcpp::as<arma::mat>(target);
  if (mat_x.n_rows != m)
    Rcpp::stop("Matrix row count in hrf_whiten_apply_cpp must match length(rows)");
  if (m >= 1 && mat_x.n_cols > 0) whitener.apply(mat_x);
  return Rcpp::wrap(mat_x);
}

// Compute U = L^T * bt = W_unit^{-T} * bt in O(m * p * M) so U^T * U = B * Sigma_sub * B^T
// on both contiguous and gapped AR(1)/AR(p) grids without forming an m x m matrix.
// [[Rcpp::export]]
arma::mat hrf_ar_operator_factor_cpp(const arma::mat& bt, const Rcpp::IntegerVector& rows,
    double rho = 0.0,
    Rcpp::Nullable<Rcpp::NumericVector> ar_phi = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> ar_acf = R_NilValue,
    double sig_u2 = 1.0) {
  if (!std::isfinite(rho) || std::abs(rho) >= 1.0)
    Rcpp::stop("Invalid AR(1) coefficient in hrf_ar_operator_factor_cpp");
  arma::uvec keep = hrf_obs_rows(rows);
  if (bt.n_rows != keep.n_elem)
    Rcpp::stop("Matrix row count in hrf_ar_operator_factor_cpp must match length(rows)");
  arma::mat out = bt;
  if (keep.n_elem <= 1 || out.n_cols == 0) return out;
  ArWhitener whitener;
  whitener.init(keep, rho, ar_phi, ar_acf, sig_u2, true);
  whitener.apply_inv_transpose(out);
  return out;
}
