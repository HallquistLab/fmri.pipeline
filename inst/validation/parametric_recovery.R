# Independent forward generator and analysis helpers for the parametric study.
# No package HRF, event-design, or convolution code enters forward generation.
# TR is fixed at one second; all amplitudes are signed integrated event mass.

pr_scenarios <- function() {
  base <- list(n = 240L, voxels = 12L, runs = 6L, beta = .6,
    eta_sd = .3, eta_rho = 0, eta_shared = .35, noise_sd = .04,
    noise_rho = .5, confounded = FALSE, spacing = c(9, 12))
  changes <- list(no_trial_variation = list(eta_sd = 0),
    moderate_trial_variation = list(), large_trial_variation = list(eta_sd = .9),
    dependent_trial_variation = list(eta_sd = .9, eta_rho = .6),
    active_null = list(beta = 0),
    confounded_null = list(beta = 0, confounded = TRUE),
    confounded_effect = list(confounded = TRUE))
  lapply(changes, function(x) utils::modifyList(base, x))
}

pr_scale <- function(x) as.numeric(scale(x))
pr_center <- function(x) x - rep(colMeans(x), each = nrow(x))
pr_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
pr_ar1 <- function(n, rho, sd) {
  # Burn-in also handles zero variance without consuming a different RNG stream.
  as.numeric(stats::filter(stats::rnorm(n + 200L) * sd * sqrt(1 - rho^2),
    rho, method = "recursive"))[-seq_len(200L)]
}
pr_hrf <- function(t, peak) {
  mass <- stats::pgamma(32, peak) - stats::pgamma(32, 16) / 6
  y <- (stats::dgamma(pmax(t, 0), peak) - stats::dgamma(pmax(t, 0), 16) / 6) / mass
  y[t < 0 | t > 32] <- 0
  y
}
pr_deposit <- function(n, onsets) {
  b <- matrix(0, n, length(onsets))
  for (e in seq_along(onsets)) {
    lo <- floor(onsets[e]); f <- onsets[e] - lo
    if (lo >= 0 && lo < n) b[lo + 1L, e] <- 1 - f
    if (lo + 1 >= 0 && lo + 1 < n) b[lo + 2L, e] <- f
  }
  b
}
pr_discrete_forward <- function(n, h) {
  lag <- outer(0:(n - 1L), 0:(n - 1L), "-")
  out <- matrix(0, n, n); valid <- lag >= 0 & lag < length(h)
  out[valid] <- h[lag[valid] + 1L]
  out
}

pr_generate <- function(config, seed) {
  stopifnot(config$n >= 180, config$voxels %% 2L == 0L, config$runs >= 3L)
  set.seed(seed)
  n <- config$n; v <- config$voxels
  ids <- paste0("v", seq_len(v)); runs <- paste0("run", seq_len(config$runs))
  groups <- rep(1:2, each = v / 2L)
  gain <- stats::runif(v, .75, 1.25) * rep(c(1, 1, 1, -1), length.out = v)
  peak <- c(5, 7)[groups]
  kernels <- vapply(peak, function(p) { h <- pr_hrf(0:32, p); h / sum(h) }, numeric(33))
  colnames(kernels) <- ids
  out <- list(config = config, gain = gain, kernels = kernels,
    groups = stats::setNames(paste0("parcel", groups), ids))
  fields <- c("bold", "events", "truth", "eta_truth", "signal", "eta_signal",
    "discrete_signal", "discrete_eta_signal", "noise", "trials", "amplitudes", "readout")
  for (field in fields) out[[field]] <- stats::setNames(vector("list", length(runs)), runs)
  for (r in seq_along(runs)) {
    feedback <- 20 + cumsum(stats::runif(n, config$spacing[1], config$spacing[2]))
    feedback <- feedback[feedback < n - 35]; m <- length(feedback)
    cue <- pr_scale(stats::rnorm(m)); rt <- pr_scale(stats::rnorm(m))
    surprise <- pr_scale(abs(stats::rnorm(m)))
    drive <- stats::rnorm(m) + if (config$confounded) .6 * cue + .5 * rt + .5 * surprise else 0
    rpe <- pr_scale(as.numeric(stats::filter(drive, .35, method = "recursive")))
    previous <- pr_scale(c(0, head(rpe, -1L)))
    events <- data.frame(onset = c(feedback - 5 - .3 * tanh(cue),
      feedback - 2.5 - .25 * tanh(rt), feedback),
      condition = rep(c("cue", "response", "feedback"), each = m))
    # Exogenous neural innovations are NOT orthogonalized against realized RPE.
    # Their accidental finite-sample association is part of inference uncertainty.
    common_eta <- pr_ar1(m, config$eta_rho, config$eta_sd)
    eta <- vapply(seq_len(v), function(j) config$eta_shared * common_eta +
      sqrt(1 - config$eta_shared^2) * pr_ar1(m, config$eta_rho, config$eta_sd), numeric(m))
    confound <- if (config$confounded) .35 * surprise + .25 * previous else rep(0, m)
    a <- rbind(matrix(.7 + .5 * cue, m, v), matrix(-.35 + .4 * rt, m, v),
      matrix(1 + config$beta * rpe + confound, m, v) + eta)
    a <- sweep(a, 2L, gain, "*")
    a_eta <- rbind(matrix(0, 2L * m, v), sweep(eta, 2L, gain, "*"))
    b <- pr_deposit(n, events$onset)
    # Normalize the linear onset readout only for fractional deposition geometry.
    # Cue/response impulses are >2 TR away and cannot enter the ideal readout.
    s <- t(pr_deposit(n, feedback)); s <- s / rowSums(s^2)
    signal <- eta_signal <- discrete <- discrete_eta <- matrix(0, n, v)
    for (g in unique(groups)) {
      j <- which(groups == g)
      continuous <- vapply(events$onset, function(at) pr_hrf(0:(n - 1L) - at, c(5, 7)[g]), numeric(n))
      discrete_design <- pr_discrete_forward(n, kernels[, j[1L]]) %*% b
      signal[, j] <- continuous %*% a[, j, drop = FALSE]
      eta_signal[, j] <- continuous %*% a_eta[, j, drop = FALSE]
      discrete[, j] <- discrete_design %*% a[, j, drop = FALSE]
      discrete_eta[, j] <- discrete_design %*% a_eta[, j, drop = FALSE]
    }
    shared_noise <- pr_ar1(n, config$noise_rho, config$noise_sd)
    noise <- vapply(seq_len(v), function(j) .35 * shared_noise + sqrt(1 - .35^2) *
      pr_ar1(n, config$noise_rho, config$noise_sd), numeric(n))
    truth <- b %*% a; eta_truth <- b %*% a_eta
    for (field in c("signal", "eta_signal", "noise", "truth", "eta_truth")) {
      value <- get(field); colnames(value) <- ids; out[[field]][[r]] <- value
    }
    colnames(discrete) <- colnames(discrete_eta) <- ids
    out$discrete_signal[[r]] <- discrete; out$discrete_eta_signal[[r]] <- discrete_eta
    out$bold[[r]] <- 100 + signal + noise
    colnames(out$bold[[r]]) <- ids
    out$events[[r]] <- events; out$amplitudes[[r]] <- a
    out$readout[[r]] <- s
    out$trials[[r]] <- data.frame(onset = feedback, rpe = rpe, cue = cue,
      rt = rt, surprise = surprise, previous = previous)
  }
  # Reuse the small regression factorizations for every method and lambda.
  out$regression_weights <- lapply(out$trials, function(trials) lapply(c(FALSE, TRUE), function(adjusted) {
    design <- if (adjusted) stats::model.matrix(~ rpe + cue + rt + surprise + previous, trials) else
      stats::model.matrix(~ rpe, trials)
    if (nrow(design) <= ncol(design) || qr(design)$rank != ncol(design)) stop("Unidentified trial regression")
    as.numeric(qr.solve(design, diag(nrow(design)))[2L, ])
  }))
  out
}

# Run-level regression followed by a t interval over run slopes. In particular,
# trials and voxels are NOT treated as independent inferential replicates.
# Cross-fitted shape selection can couple runs: coverage is audited empirically.
pr_run_slope <- function(y, trials, adjusted = TRUE) {
  design <- if (adjusted) stats::model.matrix(~ rpe + cue + rt + surprise + previous, trials) else
    stats::model.matrix(~ rpe, trials)
  if (nrow(design) <= ncol(design) || qr(design)$rank != ncol(design)) stop("Unidentified trial regression")
  as.numeric(stats::lm.fit(design, y)$coefficients[2L, ])
}
pr_inference <- function(slopes, target) {
  runs <- nrow(slopes); estimate <- colMeans(slopes)
  se <- apply(slopes, 2L, stats::sd) / sqrt(runs)
  critical <- stats::qt(.975, runs - 1L)
  statistic <- estimate / se
  statistic[se == 0 & estimate == 0] <- 0
  p <- 2 * stats::pt(-abs(statistic), runs - 1L)
  list(estimate = estimate, se = se, lower = estimate - critical * se,
    upper = estimate + critical * se, p = p,
    covered = abs(estimate - target) <= critical * se + 1e-12 * pmax(1, abs(target)),
    rejected = p < .05)
}

pr_operator <- function(n, h, rho, lambda, cache) {
  key <- paste(n, sprintf("%.17g", rho), sprintf("%.17g", lambda),
    paste(sprintf("%.17g", h), collapse = ","), sep = ":")
  if (!exists(key, cache, inherits = FALSE)) {
    d <- fmri.pipeline:::hrf_inverse_operator(n, h, seq_len(n), rho,
      matrix(0, n, 0), list(backend = "cpp", lambda = lambda, penalty = "diff2",
        ridge_floor = 1e-8, tail_padding = 0L))
    assign(key, d, cache)
  }
  get(key, cache, inherits = FALSE)
}
pr_apply <- function(y, eta_signal, h, rho, lambda, cache) {
  n <- nrow(y); v <- ncol(y)
  if (length(rho) == 1L) rho <- rep(rho, v)
  activity <- eta <- matrix(NA_real_, n, v)
  keys <- paste(apply(h, 2L, function(x) paste(sprintf("%.17g", x), collapse = ",")), rho)
  for (key in unique(keys)) {
    j <- which(keys == key)
    d <- pr_operator(n, h[, j[1L]], rho[j[1L]], lambda, cache)
    activity[, j] <- d %*% pr_center(y[, j, drop = FALSE])
    eta[, j] <- d %*% pr_center(eta_signal[, j, drop = FALSE])
  }
  list(activity = activity, eta = eta)
}

pr_cor <- function(x, y) {
  x <- pr_center(x); y <- pr_center(y)
  denominator <- sqrt(colSums(x^2) * colSums(y^2))
  value <- colSums(x * y) / denominator
  value[denominator < 1e-14] <- NA_real_; value
}
pr_run_metrics <- function(activity, eta, data, r) {
  truth <- data$truth[[r]]; interior <- 33:(nrow(truth) - 32L)
  x <- pr_center(truth[interior, , drop = FALSE]); y <- pr_center(activity[interior, , drop = FALSE])
  amplitude <- data$readout[[r]] %*% activity
  a <- tail(data$amplitudes[[r]], nrow(amplitude))
  eta_true <- data$readout[[r]] %*% data$eta_truth[[r]]
  eta_recovered <- data$readout[[r]] %*% eta
  eta_observed <- amplitude - (a - eta_true)
  eta_true <- pr_center(eta_true); eta_recovered <- pr_center(eta_recovered)
  eta_observed <- pr_center(eta_observed)
  eta_energy <- colSums(eta_true^2)
  eta_gain <- colSums(eta_true * eta_recovered) / eta_energy
  eta_gain[eta_energy < 1e-14] <- NA_real_
  eta_nrmse <- sqrt(colSums((eta_observed - eta_true)^2) / eta_energy)
  eta_nrmse[eta_energy < 1e-14] <- NA_real_
  c(neural_r = pr_mean(pr_cor(x, y)),
    neural_nrmse = mean(sqrt(colMeans((x - y)^2) / colMeans(x^2))),
    amplitude_rmse = mean(sqrt(colMeans((pr_center(amplitude) - pr_center(a))^2)) / abs(data$gain)),
    amplitude_r = pr_mean(pr_cor(amplitude, a)),
    eta_r = pr_mean(pr_cor(eta_true, eta_recovered)), eta_gain = pr_mean(eta_gain),
    eta_observed_r = pr_mean(pr_cor(eta_true, eta_observed)), eta_observed_nrmse = pr_mean(eta_nrmse),
    rpe_r = mean(sign(data$gain) * pr_cor(amplitude,
      matrix(data$trials[[r]]$rpe, nrow(amplitude), ncol(amplitude)))))
}

pr_dataset_metrics <- function(recovered, data) {
  run_metrics <- vapply(seq_along(recovered), function(r)
    pr_run_metrics(recovered[[r]]$activity, recovered[[r]]$eta, data, r), numeric(9))
  out <- rowMeans(run_metrics, na.rm = TRUE); out[!is.finite(out)] <- NA_real_
  for (adjusted in c(FALSE, TRUE)) {
    slopes <- true_slopes <- matrix(0, length(recovered), length(data$gain))
    for (r in seq_along(recovered)) {
      amplitude <- data$readout[[r]] %*% recovered[[r]]$activity
      true_amplitude <- tail(data$amplitudes[[r]], nrow(amplitude))
      weights <- data$regression_weights[[r]][[1L + as.integer(adjusted)]]
      slopes[r, ] <- crossprod(weights, amplitude)
      true_slopes[r, ] <- crossprod(weights, true_amplitude)
    }
    inference <- pr_inference(slopes, data$config$beta * data$gain)
    target <- data$config$beta * data$gain
    measures <- c(beta_bias = mean(inference$estimate / data$gain - data$config$beta),
      beta_rmse = sqrt(mean((inference$estimate / data$gain - data$config$beta)^2)),
      sample_slope_rmse = sqrt(mean(((inference$estimate - colMeans(true_slopes)) / data$gain)^2)),
      coverage = mean(inference$covered), mean_se = mean(inference$se / abs(data$gain)),
      rejection = mean(inference$rejected),
      false_positive = if (data$config$beta == 0) mean(inference$rejected) else NA_real_,
      power = if (data$config$beta != 0) mean(inference$rejected) else NA_real_,
      sign_error = if (data$config$beta != 0) mean(inference$estimate * target < 0) else NA_real_,
      wrong_sign_rejection = if (data$config$beta != 0)
        mean(inference$rejected & inference$estimate * target < 0) else NA_real_)
    out <- c(out, stats::setNames(measures, paste0(if (adjusted) "adjusted_" else "unadjusted_", names(measures))))
  }
  out
}

# Every row is one independent dataset after averaging correlated voxels/runs.
pr_summarize <- function(metrics, probability_bounds = TRUE) {
  ids <- c("scenario", "replicate", "method", "lambda")
  measures <- setdiff(names(metrics), ids)
  groups <- split(metrics, interaction(metrics$scenario, metrics$method, metrics$lambda, drop = TRUE))
  do.call(rbind, lapply(groups, function(d) do.call(rbind, lapply(measures, function(measure) {
    x <- d[[measure]]; x <- x[is.finite(x)]; n <- length(x)
    se <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
    margin <- if (n > 1L) stats::qt(.975, n - 1L) * se else NA_real_
    lower <- if (n) mean(x) - margin else NA_real_
    upper <- if (n) mean(x) + margin else NA_real_
    probability <- grepl("_(coverage|rejection|false_positive|power|sign_error|wrong_sign_rejection)$", measure)
    # Cluster fractions are not binomial counts of independent voxels. Normal MC
    # intervals use independent datasets; at an all-zero/all-one boundary use a
    # conservative one-sided bound on the probability of any non-boundary dataset.
    if (probability_bounds && probability && n > 1L) {
      lower <- max(0, lower); upper <- min(1, upper)
      if (all(x == 0)) upper <- 1 - .05^(1 / n)
      if (all(x == 1)) lower <- .05^(1 / n)
    }
    data.frame(d[1L, c("scenario", "method", "lambda")], metric = measure,
      n = n, mean = if (n) mean(x) else NA_real_, mc_se = se,
      lower = lower, upper = upper)
  }))))
}
