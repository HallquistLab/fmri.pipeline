# Validation-only extensions. Preserve the preceding study and its source hashes.
source(file.path("inst", "validation", "run_parametric_recovery.R"), local = TRUE)

cr_scenarios <- function() {
  out <- pr_scenarios()
  out$off_library_effect <- utils::modifyList(out$confounded_effect, list(peaks = c(4.3, 8.3)))
  out$off_library_null <- utils::modifyList(out$confounded_null, list(peaks = c(4.3, 8.3)))
  out
}

cr_generate <- function(config, seed) {
  data <- pr_generate(config, seed)
  if (is.null(config$peaks)) return(data)
  group <- rep(1:2, each = config$voxels / 2L)
  data$kernels <- vapply(config$peaks[group], function(p) {
    h <- pr_hrf(0:32, p); h / sum(h)
  }, numeric(33))
  colnames(data$kernels) <- colnames(data$bold[[1L]])
  # Same independently generated amplitudes, schedule and noise; change HRF only.
  for (r in seq_along(data$bold)) {
    n <- config$n; m <- nrow(data$trials[[r]])
    a <- data$amplitudes[[r]]
    eta <- rbind(matrix(0, 2L * m, config$voxels), data$readout[[r]] %*% data$eta_truth[[r]])
    b <- pr_deposit(n, data$events[[r]]$onset)
    for (g in 1:2) {
      j <- which(group == g)
      continuous <- vapply(data$events[[r]]$onset, function(at)
        pr_hrf(0:(n - 1L) - at, config$peaks[g]), numeric(n))
      discrete <- pr_discrete_forward(n, data$kernels[, j[1L]]) %*% b
      data$signal[[r]][, j] <- continuous %*% a[, j, drop = FALSE]
      data$eta_signal[[r]][, j] <- continuous %*% eta[, j, drop = FALSE]
      data$discrete_signal[[r]][, j] <- discrete %*% a[, j, drop = FALSE]
      data$discrete_eta_signal[[r]][, j] <- discrete %*% eta[, j, drop = FALSE]
    }
    data$bold[[r]] <- 100 + data$signal[[r]] + data$noise[[r]]
  }
  data
}

cr_event_model <- function(trials) {
  m <- nrow(trials)
  x <- matrix(0, 3L * m, 10L, dimnames = list(NULL, c("cue_mean", "cue_value",
    "response_mean", "response_rt", "feedback_mean", "feedback_rpe",
    "feedback_cue", "feedback_rt", "feedback_surprise", "feedback_previous")))
  x[seq_len(m), 1:2] <- cbind(1, trials$cue)
  x[m + seq_len(m), 3:4] <- cbind(1, trials$rt)
  x[2L * m + seq_len(m), 5:10] <- stats::model.matrix(~ rpe + cue + rt + surprise + previous, trials)
  x
}

cr_interpolate <- function(n, onsets, lags) {
  lapply(onsets, function(onset) {
    times <- onset + lags
    if (any(times < 0 | times > n - 1L)) stop("Readout outside the observation grid")
    t(pr_deposit(n, times))
  })
}
cr_readouts <- function(n, feedback, recovered_sources) {
  m <- length(feedback)
  b <- pr_deposit(n, feedback)
  onset <- t(b); onset <- onset / rowSums(onset * t(b))
  area <- t(vapply(cr_interpolate(n, feedback, -1:3), colSums, numeric(n)))
  area <- area / rowSums(area * t(b))
  # Euclidean matched response on fixed -3:+5 lags. The unit template is the
  # frozen inverse's response to an isolated fitted-HRF impulse, after centering.
  template <- t(vapply(seq_len(m), function(i) {
    p <- cr_interpolate(n, feedback[i], -3:5)[[1L]]
    response <- as.numeric(p %*% recovered_sources[, 2L * m + i])
    energy <- sum(response^2)
    if (!is.finite(energy) || energy < 1e-12) stop("Unidentified response template")
    as.numeric(crossprod(response, p)) / energy
  }, numeric(n)))
  list(onset = onset, area = area, template = template)
}

# Contrast-specific least squares, allowing redundant nuisance columns. Never
# silently substitute an ordinary coefficient if the target is unidentified.
cr_contrast <- function(design, target) {
  x <- design[, target]; nuisance <- design[, -target, drop = FALSE]
  norm <- sqrt(colSums(nuisance^2)); keep <- norm > 1e-12
  nuisance <- nuisance[, keep, drop = FALSE]
  rank <- 0L; residual <- x
  if (ncol(nuisance)) {
    basis <- svd(nuisance / rep(norm[keep], each = nrow(nuisance)), nv = 0)
    rank <- sum(basis$d > 1e-8 * basis$d[1L])
    if (rank) {
      q <- basis$u[, seq_len(rank), drop = FALSE]
      residual <- x - as.numeric(q %*% crossprod(q, x))
    }
  }
  energy <- sum(residual^2)
  identified <- nrow(design) - rank >= 2L && energy > 1e-10 * max(1, sum(x^2))
  list(weights = if (identified) residual / energy else rep(NA_real_, length(x)),
    identified = identified, nuisance_rank = rank,
    information_fraction = energy / max(sum(x^2), .Machine$double.eps))
}

cr_run <- function(data, fit, r, source, lambda, cache) {
  stopifnot(source %in% c("learned", "oracle", "matched_noiseless"))
  n <- nrow(data$bold[[r]]); v <- ncol(data$bold[[r]]); m <- nrow(data$trials[[r]])
  learned <- source == "learned"
  h <- if (learned) fit$kernels[, fit$selections[[r]]$kernel_index, drop = FALSE] else data$kernels
  rho <- if (learned) fit$noise[[r]]$rho else rep(data$config$noise_rho, v)
  discrete <- source == "matched_noiseless"
  y <- if (discrete) data$discrete_signal[[r]] else data$signal[[r]] + data$noise[[r]]
  activity <- matrix(NA_real_, n, v)
  out <- lapply(c("onset", "area", "template"), function(x) list(amplitude = matrix(NA_real_, m, v),
    ordinary = rep(NA_real_, v), transformed = rep(NA_real_, v),
    information = rep(NA_real_, v), nuisance_rank = rep(NA_real_, v),
    self_gain = rep(NA_real_, v), other_feedback_leakage = rep(NA_real_, v),
    cue_response_leakage = rep(NA_real_, v)))
  names(out) <- c("onset", "area", "template")
  model <- cr_event_model(data$trials[[r]])
  keys <- paste(apply(h, 2L, function(x) paste(sprintf("%.17g", x), collapse = ",")), sprintf("%.17g", rho))
  for (key in unique(keys)) {
    j <- which(keys == key)
    d <- pr_operator(n, h[, j[1L]], rho[j[1L]], lambda, cache)
    forward <- fmri.pipeline:::hrf_event_design_cpp(n, 1, 0, data$events[[r]]$onset,
      rep(0, 3L * m), seq_len(3L * m), 3L * m, h[, j[1L], drop = FALSE])[[1L]]
    # Shared BLAS application for observations and all unit-event source columns.
    applied <- d %*% pr_center(cbind(y[, j, drop = FALSE], forward))
    activity[, j] <- applied[, seq_along(j), drop = FALSE]
    g <- applied[, length(j) + seq_len(3L * m), drop = FALSE]
    readouts <- cr_readouts(n, data$trials[[r]]$onset, g)
    for (name in names(readouts)) {
      s <- readouts[[name]]
      a <- s %*% activity[, j, drop = FALSE]
      resolution <- s %*% g
      transformed <- resolution %*% model
      contrast <- cr_contrast(transformed, match("feedback_rpe", colnames(model)))
      out[[name]]$amplitude[, j] <- a
      out[[name]]$ordinary[j] <- crossprod(data$regression_weights[[r]][[2L]], a)
      out[[name]]$transformed[j] <- crossprod(contrast$weights, a)
      out[[name]]$information[j] <- contrast$information_fraction
      out[[name]]$nuisance_rank[j] <- contrast$nuisance_rank
      out[[name]]$self_gain[j] <- mean(diag(resolution[, 2L * m + seq_len(m), drop = FALSE]))
      feedback <- abs(resolution[, 2L * m + seq_len(m), drop = FALSE]); diag(feedback) <- 0
      out[[name]]$other_feedback_leakage[j] <- mean(rowSums(feedback))
      out[[name]]$cue_response_leakage[j] <- mean(rowSums(abs(resolution[, seq_len(2L * m), drop = FALSE])))
    }
  }
  list(activity = activity, readouts = out,
    deterministic = source == "matched_noiseless" && data$config$eta_sd == 0)
}

cr_metrics <- function(runs, data, readout, estimator) {
  v <- length(data$gain); nr <- length(runs)
  slopes <- do.call(rbind, lapply(runs, function(x) x$readouts[[readout]][[estimator]]))
  values <- vapply(seq_along(runs), function(r) {
    result <- runs[[r]]; z <- result$readouts[[readout]]
    truth <- tail(data$amplitudes[[r]], nrow(z$amplitude))
    interior <- 33:(nrow(data$bold[[r]]) - 32L)
    c(amplitude_r = pr_mean(pr_cor(z$amplitude, truth)),
      amplitude_rmse = mean(sqrt(colMeans((pr_center(z$amplitude) - pr_center(truth))^2)) / abs(data$gain)),
      neural_r = pr_mean(pr_cor(result$activity[interior, , drop = FALSE], data$truth[[r]][interior, , drop = FALSE])),
      self_gain = mean(z$self_gain), other_feedback_leakage = mean(z$other_feedback_leakage),
      cue_response_leakage = mean(z$cue_response_leakage),
      transformed_information_fraction = mean(z$information), transformed_nuisance_rank = mean(z$nuisance_rank))
  }, numeric(8L))
  c(cr_coefficient_metrics(slopes, data, all(vapply(runs, function(x) isTRUE(x$deterministic), logical(1)))),
    rowMeans(values))
}

cr_coefficient_metrics <- function(slopes, data, deterministic = FALSE) {
  nr <- nrow(slopes)
  # All runs must have an identified coefficient; report missingness explicitly.
  usable <- colSums(is.finite(slopes)) == nr
  inference <- if (any(usable)) pr_inference(slopes[, usable, drop = FALSE],
    data$config$beta * data$gain[usable]) else NULL
  target <- data$config$beta
  coefficient <- if (any(usable)) c(
    beta = mean(inference$estimate / data$gain[usable]),
    beta_bias = mean(inference$estimate / data$gain[usable] - target),
    beta_rmse = sqrt(mean((inference$estimate / data$gain[usable] - target)^2)),
    coefficient_coverage = mean(inference$covered),
    coefficient_false_positive = if (target == 0) mean(inference$rejected) else NA_real_,
    coefficient_power = if (target != 0) mean(inference$rejected) else NA_real_,
    coefficient_sign_error = if (target != 0) mean(inference$estimate * data$gain[usable] * target < 0) else NA_real_,
    mean_se = mean(inference$se / abs(data$gain[usable]))) else
      stats::setNames(rep(NA_real_, 8L), c("beta", "beta_bias", "beta_rmse", "coefficient_coverage",
        "coefficient_false_positive", "coefficient_power", "coefficient_sign_error", "mean_se"))
  # Deterministic matched references check numerical recovery, not sampling
  # calibration of intervals whose widths would be dominated by roundoff.
  if (deterministic)
    coefficient[c("coefficient_coverage", "coefficient_false_positive", "coefficient_power", "mean_se")] <- NA_real_
  c(coefficient, identified_fraction = mean(usable))
}

# Small benchmark only: a parametric forward GLS model produces coefficients,
# not independent trial amplitudes or an alternative neural reconstruction.
cr_forward <- function(data, fit, source) {
  slopes <- matrix(NA_real_, length(data$bold), length(data$gain))
  for (r in seq_along(data$bold)) {
    n <- nrow(data$bold[[r]]); m <- nrow(data$trials[[r]])
    h <- if (source == "learned") fit$kernels[, fit$selections[[r]]$kernel_index, drop = FALSE] else data$kernels
    rho <- if (source == "learned") fit$noise[[r]]$rho else rep(data$config$noise_rho, length(data$gain))
    y <- if (source == "matched_noiseless") data$discrete_signal[[r]] else data$signal[[r]] + data$noise[[r]]
    model <- cr_event_model(data$trials[[r]])
    keys <- paste(apply(h, 2L, function(x) paste(sprintf("%.17g", x), collapse = ",")), sprintf("%.17g", rho))
    for (key in unique(keys)) {
      j <- which(keys == key); ar <- rho[j[1L]]
      forward <- fmri.pipeline:::hrf_event_design_cpp(n, 1, 0, data$events[[r]]$onset,
        rep(0, 3L * m), seq_len(3L * m), 3L * m, h[, j[1L], drop = FALSE])[[1L]]
      x <- cbind(intercept = 1, forward %*% model)
      whiten <- function(a) rbind(sqrt(1 - ar^2) * a[1L, , drop = FALSE],
        a[-1L, , drop = FALSE] - ar * a[-n, , drop = FALSE])
      contrast <- cr_contrast(whiten(x), match("feedback_rpe", colnames(x)))
      slopes[r, j] <- crossprod(contrast$weights, whiten(y[, j, drop = FALSE]))
    }
  }
  c(cr_coefficient_metrics(slopes, data, source == "matched_noiseless" && data$config$eta_sd == 0),
    stats::setNames(rep(NA_real_, 8L), c("amplitude_r", "amplitude_rmse", "neural_r", "self_gain",
      "other_feedback_leakage", "cue_response_leakage", "transformed_information_fraction", "transformed_nuisance_rank")))
}

cr_schedule <- function() {
  learned <- expand.grid(source = "learned", lambda = c(.0001, .001, .01, .1),
    readout = c("onset", "area", "template"), estimator = c("ordinary", "transformed"), stringsAsFactors = FALSE)
  reference <- expand.grid(source = c("oracle", "matched_noiseless"), lambda = .01,
    readout = c("onset", "area", "template"), estimator = c("ordinary", "transformed"), stringsAsFactors = FALSE)
  benchmark <- data.frame(source = c("learned", "oracle", "matched_noiseless"), lambda = .01,
    readout = "forward", estimator = "gls")
  out <- rbind(learned, reference, benchmark)
  out$method <- paste(out$source, out$readout, out$estimator, sep = "/")
  out
}
