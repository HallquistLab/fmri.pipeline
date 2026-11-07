// Validation prototype; compiled by sourceCpp, not installed as a package API.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>

namespace {
struct WavePlan {
  arma::mat a, vectors, laplacian;
  std::vector<arma::mat> inverse;
  arma::uvec eigen_group;
  arma::uword p, voxels;
  bool mixture, quadratic, complete_graph;
  double rho, ridge, lambda, spatial;
};
arma::mat difference(const arma::mat& x) {
  return x.rows(1, x.n_rows - 1) - x.rows(0, x.n_rows - 2);
}
arma::mat adjoint_difference(const arma::mat& x) {
  arma::mat out(x.n_rows + 1, x.n_cols, arma::fill::zeros);
  out.rows(0, x.n_rows - 1) -= x;
  out.rows(1, x.n_rows) += x;
  return out;
}
arma::mat analysis(const arma::mat& x, const WavePlan& plan) {
  if (!plan.mixture) return difference(x);
  return arma::join_cols(x.rows(0, plan.p - 1), difference(x.rows(plan.p, 2 * plan.p - 1)));
}
arma::mat adjoint(const arma::mat& x, const WavePlan& plan) {
  if (!plan.mixture) return adjoint_difference(x);
  return arma::join_cols(x.rows(0, plan.p - 1), adjoint_difference(x.rows(plan.p, 2 * plan.p - 2)));
}
arma::mat total(const arma::mat& x, const WavePlan& plan) {
  return plan.mixture ? x.rows(0, plan.p - 1) + x.rows(plan.p, 2 * plan.p - 1) : x;
}
arma::mat apply_inverse(const arma::mat& rhs, const WavePlan& plan) {
  if (plan.spatial == 0) return plan.inverse[0] * rhs;
  if (plan.complete_graph) {
    arma::vec average = arma::mean(rhs, 1);
    arma::mat centered = rhs.each_col() - average;
    arma::mat answer = plan.inverse[1] * centered;
    answer.each_col() += plan.inverse[0] * average;
    return answer;
  }
  arma::mat rotated = rhs * plan.vectors, answer(rotated.n_rows, rotated.n_cols);
  for (arma::uword g = 0; g < plan.inverse.size(); ++g) {
    arma::uvec cols = arma::find(plan.eigen_group == g);
    answer.cols(cols) = plan.inverse[g] * rotated.cols(cols);
  }
  return answer * plan.vectors.t();
}
}

// One factorization per HRF/mask/penalty and distinct spatial eigenvalue.
// No event onsets, covariates, or neural truth are accepted by the solver.
// [[Rcpp::export]]
SEXP wave_prepare_cpp(const arma::mat& design, std::string method,
    const arma::mat& laplacian, double spatial = 0, double lambda = 0.01,
    double admm_rho = 0.01, double ridge = 1e-8) {
  if (design.n_rows < 3 || design.n_cols < 3 || !design.is_finite() ||
      laplacian.n_rows < 1 || laplacian.n_rows != laplacian.n_cols ||
      !laplacian.is_finite() || !arma::approx_equal(laplacian, laplacian.t(), "absdiff", 1e-10) ||
      !std::isfinite(spatial) || spatial < 0 || !std::isfinite(lambda) || lambda < 0 ||
      !std::isfinite(admm_rho) || admm_rho <= 0 || !std::isfinite(ridge) || ridge <= 0 ||
      (method != "diff2" && method != "tv" && method != "mixture"))
    Rcpp::stop("Invalid waveform solver design or settings");
  Rcpp::XPtr<WavePlan> plan(new WavePlan(), true);
  plan->a = design; plan->laplacian = laplacian; plan->p = design.n_cols; plan->voxels = laplacian.n_rows;
  plan->mixture = method == "mixture";
  plan->quadratic = method == "diff2"; plan->rho = admm_rho;
  plan->lambda = lambda; plan->ridge = ridge; plan->spatial = spatial;
  arma::mat complete = arma::eye<arma::mat>(laplacian.n_rows, laplacian.n_rows) -
    arma::ones<arma::mat>(laplacian.n_rows, laplacian.n_rows) / laplacian.n_rows;
  plan->complete_graph = spatial > 0 && arma::approx_equal(laplacian, complete, "absdiff", 1e-12);
  arma::vec eigen;
  if (plan->complete_graph) eigen = arma::vec({0,1});
  else if (!arma::eig_sym(eigen, plan->vectors, laplacian) || eigen.min() < -1e-9)
    Rcpp::stop("Spatial Laplacian must be positive semidefinite");
  if (spatial == 0) { eigen.zeros(); plan->vectors.eye(laplacian.n_rows, laplacian.n_cols); }
  arma::mat g = design.t() * design, base = g;
  arma::uword p = plan->p;
  if (plan->mixture) base = arma::join_cols(arma::join_rows(g, g), arma::join_rows(g, g));
  arma::mat identity = arma::eye<arma::mat>(base.n_rows, base.n_rows);
  arma::mat d;
  if (plan->quadratic) {
    d.zeros(p - 2, p);
    for (arma::uword i = 0; i < p - 2; ++i) { d(i,i) = 1; d(i,i+1) = -2; d(i,i+2) = 1; }
    base += lambda * d.t() * d;
  } else base += admm_rho * adjoint(analysis(identity, *plan), *plan);
  base.diag() += ridge;
  arma::mat ctc = arma::eye<arma::mat>(p, p);
  if (plan->mixture) ctc = arma::join_cols(arma::join_rows(ctc, ctc), arma::join_rows(ctc, ctc));
  plan->eigen_group.set_size(eigen.n_elem);
  std::vector<double> unique;
  for (arma::uword i = 0; i < eigen.n_elem; ++i) {
    double e = std::max(0.0, eigen(i));
    arma::uword j = 0;
    for (; j < unique.size(); ++j) if (std::abs(e - unique[j]) < 1e-9) break;
    if (j == unique.size()) {
      arma::mat inv;
      if (!arma::inv_sympd(inv, base + spatial * e * ctc)) Rcpp::stop("Waveform factorization failed");
      plan->inverse.push_back(inv); unique.push_back(e);
    }
    plan->eigen_group(i) = j;
  }
  return plan;
}

// Batched ADMM with exact common quadratic updates. Both primal and dual
// stopping criteria must pass for EVERY voxel. Failure is returned explicitly.
// [[Rcpp::export]]
Rcpp::List wave_solve_cpp(SEXP pointer, const arma::mat& observations,
    double lambda = 0.01, int max_iterations = 12000,
    double absolute_tolerance = 1e-6, double relative_tolerance = 1e-5) {
  Rcpp::XPtr<WavePlan> plan(pointer);
  if (observations.n_rows != plan->a.n_rows || observations.n_cols != plan->voxels ||
      !observations.is_finite() || !std::isfinite(lambda) || lambda < 0 || max_iterations < 1 ||
      !std::isfinite(absolute_tolerance) || absolute_tolerance <= 0 ||
      !std::isfinite(relative_tolerance) || relative_tolerance <= 0)
    Rcpp::stop("Invalid waveform observations or stopping settings");
  arma::mat rhs = plan->a.t() * observations;
  if (plan->mixture) rhs = arma::join_cols(rhs, rhs);
  arma::mat x = apply_inverse(rhs, *plan);
  arma::uword v = observations.n_cols;
  arma::rowvec primal(v, arma::fill::zeros), dual(v, arma::fill::zeros),
    primal_limit(v, arma::fill::zeros), dual_limit(v, arma::fill::zeros);
  arma::uvec converged(v, arma::fill::ones);
  int iterations = 0;
  if (!plan->quadratic) {
    arma::mat z = analysis(x, *plan), u(z.n_rows, v, arma::fill::zeros);
    converged.zeros();
    for (int i = 1; i <= max_iterations; ++i) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      x = apply_inverse(rhs + plan->rho * adjoint(z - u, *plan), *plan);
      arma::mat dx = analysis(x, *plan), previous = z, value = dx + u;
      z = arma::sign(value) % arma::clamp(arma::abs(value) - lambda / plan->rho, 0, arma::datum::inf);
      u += dx - z;
      primal = arma::sqrt(arma::sum(arma::square(dx - z), 0));
      dual = plan->rho * arma::sqrt(arma::sum(arma::square(adjoint(z - previous, *plan)), 0));
      arma::rowvec nx = arma::sqrt(arma::sum(arma::square(dx), 0)), nz = arma::sqrt(arma::sum(arma::square(z), 0));
      primal_limit = std::sqrt((double)z.n_rows) * absolute_tolerance + relative_tolerance * arma::max(nx, nz);
      dual_limit = std::sqrt((double)x.n_rows) * absolute_tolerance +
        relative_tolerance * plan->rho * arma::sqrt(arma::sum(arma::square(adjoint(u, *plan)), 0));
      for (arma::uword j = 0; j < v; ++j) converged(j) = primal(j) <= primal_limit(j) && dual(j) <= dual_limit(j);
      iterations = i;
      if (arma::all(converged == 1)) break;
    }
  }
  arma::mat activity = total(x, *plan), residual = plan->a * activity - observations;
  double regularization = plan->quadratic ?
    .5 * plan->lambda * arma::accu(arma::square(difference(difference(x)))) : lambda * arma::accu(arma::abs(analysis(x, *plan)));
  double objective = .5 * arma::accu(arma::square(residual)) + regularization + .5 * plan->ridge * arma::accu(arma::square(x));
  objective += .5 * plan->spatial * arma::accu(activity % (activity * plan->laplacian));
  return Rcpp::List::create(Rcpp::Named("activity") = activity,
    Rcpp::Named("components") = x, Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged, Rcpp::Named("primal") = primal,
    Rcpp::Named("dual") = dual, Rcpp::Named("primal_limit") = primal_limit,
    Rcpp::Named("dual_limit") = dual_limit,
    Rcpp::Named("objective") = objective);
}
