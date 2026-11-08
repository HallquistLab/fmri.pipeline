#' Diagnose amplitude attenuation and event leakage for a frozen inverse
#'
#' @param deconvolution Result of deconvolve_crossrun_hrf().
#' @param events Optional run-named event tables with onset and optional condition
#'   and duration. Defaults to the HRF training event records. Every record is a
#'   separate unit-amplitude source, regardless of its condition label.
#' @param lags Readout times relative to each event onset, in seconds. Default 0.
#' @param weights Linear readout weights, default one per lag. For example,
#'   lags = -1:3 produces a five-sample window sum at TR = 1.
#' @param source_kernel_grid Optional normalized HRF candidates used to generate
#'   synthetic BOLD. Default uses each stratum's fitted HRF. Supplying alternatives
#'   evaluates sensitivity to HRF mismatch; it does not estimate the true HRF.
#' @param return_matrices Return event-by-event recovered, ideal, and excess
#'   resolution matrices, default FALSE.
#' @param max_events Explicit limit per run on the quadratic diagnostic matrix,
#'   default 2000. Use a smaller prespecified subset or explicitly raise the limit.
#' @return A list with per-event diagnostics, stratum summaries, voxel-to-stratum
#'   mappings, optional resolution matrices, and the conditional assumptions.
#' @details For source events B, convolution H, frozen inverse D, and readout S,
#'   recovered = S D H B and ideal = S B. Gain is their diagonal ratio. Excess
#'   leakage is recovered minus ideal, normalized by each receiving event's ideal
#'   diagonal. Subtracting ideal retains natural overlap of neural readout windows.
#'   Previous/next refer to chronological source records, including other event
#'   types. Zero ideal diagonals and readouts outside the acquisition grid are NA.
#'
#'   The calculation reuses one operator per fitted HRF/noise stratum and uses
#'   BLAS across source events, not a simulation or fit per voxel. It includes the
#'   actual masks, nuisance projection, padding, centering, and inverse penalty.
#'   Scaling divisors cancel in the gain ratios and are omitted. Their estimated
#'   values, noise parameters, selected HRFs, and lambda are held fixed: these are
#'   conditional linear diagnostics, not uncertainty intervals or a claim about
#'   physiological neural units. Source events of positive duration have unit
#'   height (fractional acquisition-bin occupancy), not unit integrated area.
#'   No automatic correction is applied; a scalar cannot undo a full matrix of
#'   attenuation and cross-talk. Wiener diagnostics are not supported here.
#' @export
diagnose_crossrun_recovery <- function(deconvolution, events = NULL, lags = 0,
    weights = NULL, source_kernel_grid = NULL, return_matrices = FALSE, max_events = 2000L) {
  fit <- deconvolution$hrf_fit; settings <- deconvolution$settings
  if (!inherits(fit, "crossrun_hrf") || is.null(deconvolution$scaling) ||
      is.null(settings$nuisance_mode)) stop("Use a current deconvolve_crossrun_hrf result")
  if (settings$backend == "wiener") stop("Recovery diagnostics require the cpp or r GLS inverse")
  checkmate::assert_numeric(lags, min.len = 1L, any.missing = FALSE, finite = TRUE)
  if (is.null(weights)) weights <- rep(1, length(lags))
  checkmate::assert_numeric(weights, len = length(lags), any.missing = FALSE, finite = TRUE)
  checkmate::assert_flag(return_matrices)
  checkmate::assert_integerish(max_events, lower = 1L, len = 1L)
  if (is.null(events)) events <- fit$specification$events
  if (!is.list(events) || !identical(names(events), fit$run_ids))
    stop("events must be a list in fitted run order")
  sources <- if (is.null(source_kernel_grid)) NULL else
    do.call(cbind, build_decon_kernel_grid(c(0, 1), source_kernel_grid))
  if (!is.null(sources) && !is.null(attr(source_kernel_grid, "tr")) &&
      !isTRUE(all.equal(attr(source_kernel_grid, "tr"), fit$settings$tr)))
    stop("source_kernel_grid sampling TR does not match the fit")
  all_events <- summaries <- mappings <- matrices <- list()
  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  for (r in seq_along(fit$run_ids)) {
    n <- fit$n_time[r]; ev <- events[[r]]; tr <- fit$settings$tr; offset <- fit$settings$time_offset[r]
    checkmate::assert_data_frame(ev, min.rows = 1L)
    checkmate::assert_numeric(ev$onset, any.missing = FALSE, finite = TRUE)
    if (nrow(ev) > max_events) stop("Diagnostic event matrix exceeds max_events")
    duration <- if (is.null(ev$duration)) rep(0, nrow(ev)) else ev$duration
    checkmate::assert_numeric(duration, lower = 0, len = nrow(ev), any.missing = FALSE, finite = TRUE)
    m <- nrow(ev); order <- order(ev$onset, seq_len(m)); position <- match(seq_len(m), order)
    previous <- ifelse(position > 1L, order[pmax(1L, position - 1L)], NA_integer_)
    next_event <- ifelse(position < m, order[pmin(m, position + 1L)], NA_integer_)
    readout <- matrix(0, m, n); observed <- rep(TRUE, m)
    for (j in seq_along(lags)) {
      at <- (ev$onset + lags[j] - offset) / tr
      lo <- floor(at); fraction <- at - lo
      ok <- at >= 0 & at <= n - 1L
      if (weights[j] == 0) next
      observed <- observed & ok
      index <- which(ok)
      readout[cbind(index, lo[index] + 1L)] <- readout[cbind(index, lo[index] + 1L)] + weights[j] * (1 - fraction[index])
      index <- which(ok & fraction > 0 & lo + 2L <= n)
      readout[cbind(index, lo[index] + 2L)] <- readout[cbind(index, lo[index] + 2L)] + weights[j] * fraction[index]
    }
    readout[!observed, ] <- NA_real_
    neural <- hrf_event_design_cpp(n, tr, offset, ev$onset, duration,
      seq_len(m), m, matrix(c(1, 0), ncol = 1L))[[1L]]
    ideal <- readout %*% neural
    ideal_diagonal <- diag(ideal)
    mask <- if (settings$observation_mask == "fit") fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    if (length(rows) < 3L) next
    z <- if (settings$nuisance_mode == "joint") cbind(1, seq(-1, 1, length.out = n),
      fit$specification$nuisance[[r]]) else matrix(numeric(), n, 0L)
    selected <- fit$selections[[r]]$kernel_index
    rho <- if (is.null(fit$noise)) rep(0, length(selected)) else fit$noise[[r]]$rho
    keys <- paste(selected, sprintf("%.17g", rho), sep = ":")
    valid <- deconvolution$scaling[[r]]$status == "ok"
    mappings[[r]] <- data.frame(run = fit$run_ids[r], voxel = fit$voxel_ids,
      stratum = keys, usable = valid)
    for (key in unique(keys[valid])) {
      voxels <- which(keys == key & valid); k <- selected[voxels[1L]]; ar <- rho[voxels[1L]]
      operator <- hrf_inverse_operator(n, fit$kernels[, k], rows, ar, z, settings)
      readout_operator <- readout %*% operator
      source <- if (is.null(sources)) fit$kernels[, k, drop = FALSE] else sources
      designs <- hrf_event_design_cpp(n, tr, offset, ev$onset, duration, seq_len(m), m, source)
      for (s in seq_len(ncol(source))) {
        y <- designs[[s]][rows, , drop = FALSE]
        if (settings$scaling != "none") y <- sweep(y, 2L, colMeans(y), "-")
        recovered <- readout_operator %*% y
        excess <- recovered - ideal
        divisor <- ideal_diagonal
        divisor[!is.finite(divisor) | abs(divisor) < 1e-10] <- NA_real_
        normalized <- sweep(excess, 1L, divisor, "/")
        gain <- diag(recovered) / divisor
        off_diagonal <- normalized; diag(off_diagonal) <- NA_real_
        meta <- data.frame(run = fit$run_ids[r], stratum = key, fitted_kernel = colnames(fit$kernels)[k],
          source_kernel = if (is.null(sources)) "fitted" else paste0("source_", s), noise_rho = ar,
          n_voxels = length(voxels))
        all_events[[length(all_events) + 1L]] <- cbind(meta[rep(1L, m), ],
          data.frame(event = seq_len(m), onset = ev$onset,
            condition = if (is.null(ev$condition)) NA_character_ else as.character(ev$condition),
            ideal_diagonal = ideal_diagonal, recovered_diagonal = diag(recovered), amplitude_gain = gain,
            previous_event = previous, previous_excess_leakage = normalized[cbind(seq_len(m), previous)],
            next_event = next_event, next_excess_leakage = normalized[cbind(seq_len(m), next_event)],
            mean_abs_excess_leakage = apply(abs(off_diagonal), 1L, finite_mean)))
        summaries[[length(summaries) + 1L]] <- cbind(meta,
          data.frame(n_events = m, n_valid_gains = sum(is.finite(gain)), mean_gain = finite_mean(gain),
            mean_abs_gain_error = finite_mean(abs(gain - 1)),
            mean_abs_excess_leakage = finite_mean(abs(off_diagonal))))
        if (return_matrices) matrices[[paste(fit$run_ids[r], key, s, sep = "/")]] <-
          list(recovered = recovered, ideal = ideal, excess = excess, normalized_excess = normalized)
      }
    }
  }
  list(events = if (length(all_events)) do.call(rbind, all_events) else data.frame(),
    summary = if (length(summaries)) do.call(rbind, summaries) else data.frame(),
    voxels = if (length(mappings)) do.call(rbind, mappings) else data.frame(),
    matrices = if (return_matrices) matrices else NULL,
    settings = list(lags = lags, weights = weights, inverse = settings,
      conditional_on = "HRF, covariance, lambda, masks, nuisance design, and scaling divisors",
      reference = "unit-height event sources; unit-sum HRFs; ratios within analysis units",
      automatic_correction = FALSE))
}

build_epoch_readout_matrix <- function(onsets, lag, n_time, tr, offset,
    readout_kernel = "linear", readout_width = NULL) {
  m <- length(onsets)
  if (is.null(readout_width) || length(readout_width) <= 1L) {
    quad <- epoch_kernel_quadrature(readout_kernel, readout_width, tr)
    width_scale <- rep(1, m)
  } else {
    if (length(readout_width) != m) {
      stop("readout_width must be NULL, a scalar, or a vector of length(onsets)")
    }
    checkmate::assert_numeric(readout_width, lower = .Machine$double.eps, finite = TRUE, any.missing = FALSE)
    quad <- epoch_kernel_quadrature(readout_kernel, 1.0, tr)
    width_scale <- as.numeric(readout_width)
  }
  readout <- matrix(0, m, n_time)
  weight_sum <- rep(0, m)
  for (q in seq_along(quad$offsets)) {
    at <- (onsets + lag + quad$offsets[q] * width_scale - offset) / tr
    lo <- floor(at)
    fraction <- at - lo
    ok <- is.finite(at) & at >= 0 & at <= n_time - 1L
    w_q <- quad$weights[q]
    if (!any(ok) || w_q == 0) next
    idx <- which(ok)
    readout[cbind(idx, lo[idx] + 1L)] <- readout[cbind(idx, lo[idx] + 1L)] +
      w_q * (1 - fraction[idx])
    idx2 <- which(ok & fraction > 0 & lo + 2L <= n_time)
    if (length(idx2)) {
      readout[cbind(idx2, lo[idx2] + 2L)] <- readout[cbind(idx2, lo[idx2] + 2L)] +
        w_q * fraction[idx2]
    }
    weight_sum[idx] <- weight_sum[idx] + w_q
  }
  valid <- weight_sum > 0
  if (any(valid)) {
    readout[valid, ] <- sweep(readout[valid, , drop = FALSE], 1L, weight_sum[valid], "/")
  }
  readout[!valid, ] <- NA_real_
  readout
}

solve_duration_demix_operator <- function(S_tau, operator, G_tau_unc, G_tau_cent,
    scaling_mode, demix_ridge, gain_floor) {
  m_ok <- nrow(S_tau)
  n_obs <- ncol(operator)
  K_norm <- S_tau %*% (operator %*% G_tau_cent)
  K_unc <- S_tau %*% (operator %*% G_tau_unc)
  g_unc <- diag(K_unc)
  g_safe <- sign(ifelse(g_unc == 0, 1, g_unc)) * pmax(abs(g_unc), gain_floor)
  R_unc <- sweep(K_unc, 1L, g_safe, "/")
  Q <- diag(m_ok) - matrix(1 / m_ok, m_ok, m_ok)
  lhs <- crossprod(R_unc) + demix_ridge * Q + 1e-8 * diag(m_ok)
  W_unc <- solve(lhs, sweep(t(R_unc), 2L, g_safe, "/"))
  if (scaling_mode != "none") {
    c0 <- as.numeric(S_tau %*% (operator %*% rep(1, n_obs)))
    mu_G <- colMeans(G_tau_unc)
    u0 <- as.numeric(W_unc %*% c0)
    denom <- 1 - sum(mu_G * u0)
    if (!is.finite(denom) || abs(denom) < 0.05) {
      denom <- sign(ifelse(!is.finite(denom) || denom == 0, 1, denom)) * 0.05
    }
    W_full <- W_unc + outer(u0, as.numeric(crossprod(mu_G, W_unc))) / denom
  } else {
    W_full <- W_unc
  }
  eff_gain <- if (m_ok > 1L) {
    sum(diag(Q %*% W_full %*% K_norm %*% Q)) / (m_ok - 1L)
  } else {
    1 / (1 + demix_ridge)
  }
  if (!is.finite(eff_gain) || abs(eff_gain) < gain_floor) {
    eff_gain <- 1 / (1 + demix_ridge)
  }
  eff_gain_vec <- if (m_ok > 1L) {
    g_d <- diag(Q %*% W_full %*% K_norm %*% Q) / diag(Q)
    g_d[!is.finite(g_d)] <- eff_gain
    sign(ifelse(g_d == 0, 1, g_d)) * pmax(abs(g_d), gain_floor)
  } else {
    rep(eff_gain, m_ok)
  }
  list(W = W_full, K_norm = K_norm, eff_gain = eff_gain, eff_gain_vec = eff_gain_vec)
}

compute_ar1_operator_cov <- function(B, rows, rho, center_obs = TRUE) {
  m <- ncol(B)
  M <- nrow(B)
  if (m == 0L || M == 0L) return(matrix(0, M, M))
  B_eff <- if (isTRUE(center_obs)) {
    sweep(B, 1L, rowMeans(B), "-")
  } else {
    B
  }
  if (!is.finite(rho) || abs(rho) < 1e-10 || m == 1L) {
    return(tcrossprod(B_eff))
  }
  Bt <- t(B_eff)
  V <- matrix(0, m, M)
  V[m, ] <- Bt[m, ]
  a_vec <- c(0, rho^(rows[-1L] - rows[-m]))
  for (i in (m - 1L):1L) {
    V[i, ] <- Bt[i, ] + a_vec[i + 1L] * V[i + 1L, ]
  }
  s_vec <- sqrt(pmax(0, 1 - a_vec * a_vec))
  U <- sweep(V, 1L, s_vec, "*")
  crossprod(U)
}

build_other_demix_operator <- function(S_all_0, DG_all_0, ev, o_idx, t_idx, demix_ridge) {
  ok_all <- rowSums(!is.finite(S_all_0)) == 0L
  if (!any(ok_all) || !length(o_idx)) return(NULL)
  m_all <- nrow(ev)
  cond_all <- if (is.null(ev$condition)) rep("event", m_all) else as.character(ev$condition)
  cond_levels <- unique(cond_all[ok_all])
  M_cond <- vapply(cond_levels, function(cl) as.numeric(cond_all == cl), numeric(m_all))
  if (is.null(dim(M_cond))) M_cond <- matrix(M_cond, ncol = 1L)
  S_ok <- S_all_0[ok_all, , drop = FALSE]
  R_all_0 <- S_ok %*% DG_all_0
  R_cond <- R_all_0 %*% M_cond
  lhs_cond <- crossprod(R_cond) + 1e-4 * mean(diag(crossprod(R_cond))) * diag(ncol(R_cond))
  A_mean_op <- M_cond %*% solve(lhs_cond, crossprod(R_cond, S_ok))

  ok_o <- o_idx[o_idx %in% which(ok_all)]
  if (!length(ok_o)) {
    return(A_mean_op[o_idx, , drop = FALSE])
  }
  S_o <- S_all_0[ok_o, , drop = FALSE]
  R_oo <- S_o %*% DG_all_0[, o_idx, drop = FALSE]
  R_ot <- S_o %*% DG_all_0[, t_idx, drop = FALSE]
  Res_o_op <- S_o - (S_o %*% DG_all_0) %*% A_mean_op
  if (ncol(R_ot) > 0L) {
    diag_ot <- max(mean(diag(crossprod(R_ot))), 1e-6)
    P_t <- R_ot %*% solve(crossprod(R_ot) + 0.05 * diag_ot * diag(ncol(R_ot)), t(R_ot))
    Res_o_op <- Res_o_op - P_t %*% Res_o_op
  }
  diag_oo <- max(mean(diag(crossprod(R_oo))), 1e-6)
  lhs_oo <- crossprod(R_oo) + max(demix_ridge, 0.20) * diag_oo * diag(ncol(R_oo))
  Delta_o_op <- solve(lhs_oo, crossprod(R_oo, Res_o_op))
  A_mean_op[o_idx, , drop = FALSE] + Delta_o_op
}

compute_zca_whitening_matrix <- function(Sigma_tot, eig_floor = 0.05) {
  M <- nrow(Sigma_tot)
  if (M <= 1L) return(diag(M))
  d_tot <- sqrt(pmax(diag(Sigma_tot), .Machine$double.eps))
  R_tot <- sweep(sweep(Sigma_tot, 1L, d_tot, "/"), 2L, d_tot, "/")
  R_tot <- 0.5 * (R_tot + t(R_tot))
  eig <- eigen(R_tot, symmetric = TRUE)
  inv_sqrt_vals <- 1 / sqrt(pmax(eig$values, eig_floor))
  R_inv_half <- eig$vectors %*% (inv_sqrt_vals * t(eig$vectors))
  Phi <- sweep(sweep(R_inv_half, 1L, d_tot, "*"), 2L, d_tot, "/")
  tr_mean <- mean(diag(Phi))
  if (is.finite(tr_mean) && abs(tr_mean) > 1e-6) {
    Phi <- Phi / tr_mean
  }
  Phi
}

#' Operator-based epoch de-mixing and gain-corrected trial regression across lags
#'
#' Estimate trial-by-trial correlations and regression slopes across epoch lags
#' relative to a target event while removing preceding/neighboring event leakage
#' and correcting regularization-induced amplitude attenuation.
#'
#' @param deconvolution Result of \code{deconvolve_crossrun_hrf()}.
#' @param covariates Numeric vector, matrix, or data.frame of trial-level
#'   predictors for target trials (concatenated across runs or supplied as a
#'   run-named list matching \code{deconvolution$hrf_fit$run_ids}).
#' @param events Optional run-named list of event data.frames with \code{onset}
#'   and optional \code{condition} and \code{duration}. Defaults to the HRF fit's
#'   event tables.
#' @param target_condition Optional character scalar naming the target event
#'   condition in \code{events[[r]]$condition} (e.g., \code{"feedback"}). When
#'   supplied, non-target events in \code{events[[r]]} are treated as interfering
#'   events whose leakage into the target epoch is removed when
#'   \code{demix_other_events = TRUE}.
#' @param lags Numeric vector of epoch readout lags (in seconds) relative to each
#'   target event onset.
#' @param readout_kernel Epoch readout kernel: \code{"linear"}, \code{"boxcar"},
#'   \code{"triangle"}, or \code{"gaussian"}. Default: \code{"triangle"}.
#' @param readout_width Optional kernel full width in seconds for non-linear
#'   readouts (defaults to \code{1.5 * tr}).
#' @param demix_other_events Logical; whether to estimate non-target event
#'   amplitudes at their own onsets and subtract their operator-predicted leakage
#'   into each target epoch lag. Default: \code{TRUE}.
#' @param demix_target_trials Logical; whether to use the stratum resolution
#'   matrix at each lag to de-mix adjacent target-trial overlap and correct
#'   amplitude attenuation. Default: \code{TRUE}.
#' @param demix_ridge Positive ridge stabilization parameter for unit-diagonal
#'   cross-talk inversion. Default: \code{0.05}.
#' @param gain_floor Minimum absolute operator gain used when correcting
#'   amplitude slopes. Default: \code{0.05}.
#' @param restore_scale Logical; if \code{TRUE}, multiplies each run and voxel by
#'   its recorded \code{deconvolution$scaling[[r]]$scale} divisor so that
#'   \code{beta_raw} and \code{beta_corrected} are reported in unscaled input
#'   units rather than standardized (\code{zscore}/\code{noise_sd}) units.
#'   Default: \code{FALSE}.
#' @param source_model Forward neural source model at each target lag:
#'   \code{"impulse"} (default unit impulse at \code{onset + lag}),
#'   \code{"window"} (local sustained boxcar of width \code{readout_width}
#'   centered at \code{onset + lag}), or \code{"duration"} (trial-varying
#'   sustained boxcar over \code{[onset + lag, onset + lag + d_i]} using each
#'   target trial's duration \code{d_i} from \code{trial_durations} or
#'   \code{events[[r]]$duration}, with duration-matched readout and unpenalized
#'   baseline-rate de-confounding).
#' @param trial_durations Optional numeric vector (concatenated across runs) or
#'   run-named list of non-negative target trial durations (in seconds), such as
#'   per-trial reaction times. Overrides \code{events[[r]]$duration} for target
#'   trials and defaults \code{source_model} to \code{"duration"} when
#'   \code{source_model} is not explicitly specified.
#' @param duration_scale Scaling convention for positive-duration target trials
#'   when \code{source_model = "duration"}: \code{"rate"} (default; trial
#'   amplitude represents neural firing rate / boxcar height per TR bin, so
#'   integrated neural activity scales with duration and constant baseline firing
#'   rates are de-confounded from duration variability) or \code{"integral"}
#'   (each trial boxcar is normalized to unit time-integrated neural mass across
#'   its duration).
#' @param weighting Trial regression weighting scheme across target trials and
#'   runs: \code{"precision"} (default; Empirical-Bayes / floor-stabilized
#'   inverse error-variance weights propagating per-trial operator standard
#'   errors and run noise scale), \code{"gls"} (within-run ZCA covariance
#'   whitening plus precision weighting accounting for off-diagonal trial-to-trial
#'   error correlations), or \code{"none"} (unweighted OLS). All three slope
#'   estimates (\code{beta_unweighted}, \code{beta_weighted}, \code{beta_gls})
#'   and their standard errors are always included in \code{coefficients}.
#' @param whiten_panel Logical; if \code{TRUE} and \code{return_panel = TRUE},
#'   adds within-run variance-preserving ZCA-whitened outcome
#'   (\code{decon_demixed_whitened}) and covariate (\code{<var>_whitened})
#'   columns to \code{panel} so downstream \code{mixed_by(..., weights = "precision_weight")}
#'   models account for both heteroscedastic trial precision and off-diagonal
#'   within-run trial error correlations. Default: \code{FALSE}.
#' @param return_epochs Logical; whether to include raw and de-mixed trial-by-voxel
#'   epoch matrices per lag. Default: \code{FALSE}.
#' @param return_panel Logical; whether to return a tidy MEDUSA-compatible
#'   \code{data.table} (\code{panel}) with columns \code{run}, \code{trial},
#'   \code{event_index}, \code{onset}, \code{evt_time}, \code{voxel},
#'   \code{decon_raw}, \code{decon_demixed}, \code{se_raw}, \code{se_demixed},
#'   \code{precision_weight}, \code{precision_weight_op},
#'   \code{precision_weight_raw}, \code{error_ac1}, and the covariate columns for
#'   downstream \code{mixed_by()} modeling. Default: \code{FALSE}.
#' @return A list with \code{coefficients} (per voxel, lag, and predictor term),
#'   \code{summary} (averaged across usable voxels per lag and term), optional
#'   \code{epochs}, optional \code{panel}, and \code{settings}.
#' @export
decon_epoch_regression <- function(
  deconvolution,
  covariates,
  events = NULL,
  target_condition = NULL,
  lags = seq(-2, 4, by = 1),
  readout_kernel = c("triangle", "linear", "boxcar", "gaussian"),
  readout_width = NULL,
  demix_other_events = TRUE,
  demix_target_trials = TRUE,
  demix_ridge = 0.05,
  gain_floor = 0.05,
  restore_scale = FALSE,
  source_model = c("impulse", "window", "duration"),
  trial_durations = NULL,
  duration_scale = c("rate", "integral"),
  weighting = c("precision", "gls", "none"),
  whiten_panel = FALSE,
  return_epochs = FALSE,
  return_panel = FALSE
) {
  fit <- deconvolution$hrf_fit
  settings <- deconvolution$settings
  if (!inherits(fit, "crossrun_hrf") || is.null(deconvolution$scaling) ||
      is.null(settings$nuisance_mode)) {
    stop("Use a current deconvolve_crossrun_hrf result")
  }
  if (settings$backend == "wiener") {
    stop("decon_epoch_regression requires the cpp or r GLS inverse")
  }
  source_model_supplied <- !missing(source_model)
  readout_kernel <- match.arg(readout_kernel)
  source_model <- match.arg(source_model)
  duration_scale <- match.arg(duration_scale)
  weighting <- match.arg(weighting)
  if (!is.null(trial_durations) && !source_model_supplied) {
    source_model <- "duration"
  }
  checkmate::assert_numeric(lags, min.len = 1L, any.missing = FALSE, finite = TRUE)
  checkmate::assert_string(target_condition, null.ok = TRUE)
  checkmate::assert_flag(demix_other_events)
  checkmate::assert_flag(demix_target_trials)
  checkmate::assert_number(demix_ridge, lower = 1e-6, finite = TRUE)
  checkmate::assert_number(gain_floor, lower = 1e-4, finite = TRUE)
  checkmate::assert_flag(restore_scale)
  checkmate::assert_flag(whiten_panel)
  checkmate::assert_flag(return_epochs)
  checkmate::assert_flag(return_panel)

  if (is.null(events)) events <- fit$specification$events
  if (!is.list(events) || !identical(names(events), fit$run_ids)) {
    stop("events must be a list in fitted run order")
  }

  n_runs <- length(fit$run_ids)
  target_idx_list <- vector("list", n_runs)
  other_idx_list <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    ev <- events[[r]]
    checkmate::assert_data_frame(ev, min.rows = 1L)
    checkmate::assert_numeric(ev$onset, any.missing = FALSE, finite = TRUE)
    if (is.null(target_condition)) {
      target_idx_list[[r]] <- seq_len(nrow(ev))
      other_idx_list[[r]] <- integer(0L)
    } else {
      if (is.null(ev$condition)) {
        stop("target_condition supplied, but events[[", r, "]] has no 'condition' column")
      }
      cond <- as.character(ev$condition)
      target_idx_list[[r]] <- which(!is.na(cond) & cond == target_condition)
      other_idx_list[[r]] <- which(is.na(cond) | cond != target_condition)
      if (!length(target_idx_list[[r]])) {
        stop("No events match target_condition = '", target_condition, "' in run ", fit$run_ids[r])
      }
    }
  }

  target_counts <- vapply(target_idx_list, length, integer(1L))
  cov_list <- vector("list", n_runs)
  if (is.list(covariates) && !is.data.frame(covariates)) {
    if (length(covariates) != n_runs) {
      stop("When covariates is a list, it must have one element per run")
    }
    for (r in seq_len(n_runs)) {
       x_r <- as.matrix(covariates[[r]])
      if (nrow(x_r) != target_counts[r]) {
        stop("covariates[[", r, "]] row count (", nrow(x_r),
          ") does not match target event count (", target_counts[r], ")")
      }
      if (is.null(colnames(x_r))) colnames(x_r) <- paste0("x", seq_len(ncol(x_r)))
      cov_list[[r]] <- x_r
    }
  } else {
    x_all <- as.matrix(covariates)
    if (nrow(x_all) != sum(target_counts)) {
      stop("Total rows in covariates (", nrow(x_all),
        ") must equal total target events across runs (", sum(target_counts), ")")
    }
    if (is.null(colnames(x_all))) colnames(x_all) <- paste0("x", seq_len(ncol(x_all)))
    ends <- cumsum(target_counts)
    starts <- ends - target_counts + 1L
    for (r in seq_len(n_runs)) {
      cov_list[[r]] <- x_all[starts[r]:ends[r], , drop = FALSE]
    }
  }
  term_names <- colnames(cov_list[[1L]])
  n_terms <- length(term_names)
  n_voxels <- length(fit$voxel_ids)
  n_lags <- length(lags)
  tr <- fit$settings$tr

  ends <- cumsum(target_counts)
  starts <- ends - target_counts + 1L
  dur_target_list <- vector("list", n_runs)
  if (!is.null(trial_durations)) {
    if (is.list(trial_durations) && !is.data.frame(trial_durations)) {
      if (length(trial_durations) != n_runs) {
        stop("When trial_durations is a list, it must have one element per run")
      }
      for (r in seq_len(n_runs)) {
        d_r <- as.numeric(trial_durations[[r]])
        if (length(d_r) != target_counts[r]) {
          stop("trial_durations[[", r, "]] length (", length(d_r),
            ") does not match target event count (", target_counts[r], ")")
        }
        checkmate::assert_numeric(d_r, lower = 0, any.missing = FALSE, finite = TRUE)
        dur_target_list[[r]] <- d_r
      }
    } else {
      d_all <- as.numeric(trial_durations)
      if (length(d_all) != sum(target_counts)) {
        stop("Total length of trial_durations (", length(d_all),
          ") must equal total target events across runs (", sum(target_counts), ")")
      }
      checkmate::assert_numeric(d_all, lower = 0, any.missing = FALSE, finite = TRUE)
      for (r in seq_len(n_runs)) {
        dur_target_list[[r]] <- d_all[starts[r]:ends[r]]
      }
    }
  }

  total_trials <- sum(target_counts)
  raw_epochs <- demixed_epochs <- se_raw_epochs <- se_demixed_epochs <- vector("list", n_lags)
  pw_epochs <- pw_op_epochs <- pw_raw_epochs <- error_ac1_epochs <- vector("list", n_lags)
  whitened_y_epochs <- if (isTRUE(whiten_panel)) vector("list", n_lags) else NULL
  whitened_x_epochs <- if (isTRUE(whiten_panel)) vector("list", n_lags) else NULL
  V_demix_blocks <- vector("list", n_lags)
  lag_names <- sprintf("lag_%g", lags)
  names(raw_epochs) <- names(demixed_epochs) <- names(se_raw_epochs) <- names(se_demixed_epochs) <- lag_names
  names(pw_epochs) <- names(pw_op_epochs) <- names(pw_raw_epochs) <- names(error_ac1_epochs) <- lag_names
  names(V_demix_blocks) <- lag_names
  for (li in seq_len(n_lags)) {
    na_mat <- matrix(NA_real_, total_trials, n_voxels, dimnames = list(NULL, fit$voxel_ids))
    raw_epochs[[li]] <- na_mat
    demixed_epochs[[li]] <- na_mat
    se_raw_epochs[[li]] <- na_mat
    se_demixed_epochs[[li]] <- na_mat
    pw_epochs[[li]] <- na_mat
    pw_op_epochs[[li]] <- na_mat
    pw_raw_epochs[[li]] <- na_mat
    error_ac1_epochs[[li]] <- na_mat
    if (isTRUE(whiten_panel)) {
      whitened_y_epochs[[li]] <- na_mat
      whitened_x_epochs[[li]] <- stats::setNames(
         replicate(n_terms, na_mat, simplify = FALSE),
        term_names
      )
    }
    V_demix_blocks[[li]] <- stats::setNames(vector("list", n_runs), fit$run_ids)
    for (r in seq_len(n_runs)) V_demix_blocks[[li]][[r]] <- list()
  }
  gain_by_lag_voxel <- matrix(NA_real_, n_lags, n_voxels, dimnames = list(NULL, fit$voxel_ids))
  x_centered_list <- lapply(cov_list, function(x_r) sweep(x_r, 2L, colMeans(x_r, na.rm = TRUE), "-"))
  X_raw_pool <- do.call(rbind, cov_list)

  gain_accum <- matrix(0, n_lags, n_voxels)
  gain_count <- matrix(0L, n_lags, n_voxels)
  stratum_map <- vector("list", n_runs)
  sig_eff_mat <- matrix(NA_real_, n_runs, n_voxels, dimnames = list(fit$run_ids, fit$voxel_ids))
  center_obs <- !identical(settings$scaling, "none")

  for (r in seq_len(n_runs)) {
    n <- fit$n_time[r]
    ev <- events[[r]]
    offset <- fit$settings$time_offset[r]
    t_idx <- target_idx_list[[r]]
    o_idx <- other_idx_list[[r]]
    m_all <- nrow(ev)
    m_target <- length(t_idx)
    row_span <- starts[r]:ends[r]
    duration <- if (is.null(ev$duration)) rep(0, m_all) else as.numeric(ev$duration)
    checkmate::assert_numeric(duration, lower = 0, len = m_all, any.missing = FALSE, finite = TRUE)
    if (!is.null(trial_durations)) {
      duration[t_idx] <- dur_target_list[[r]]
    }
    dur_target <- duration[t_idx]
    w_src <- if (is.null(readout_width)) 1.5 * tr else readout_width[1L]

    mask <- if (settings$observation_mask == "fit") fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    if (length(rows) < 3L) next
    z <- if (settings$nuisance_mode == "joint") {
      cbind(1, seq(-1, 1, length.out = n), fit$specification$nuisance[[r]])
    } else matrix(numeric(), n, 0L)

    selected <- fit$selections[[r]]$kernel_index
    rho <- if (is.null(fit$noise)) rep(0, length(selected)) else fit$noise[[r]]$rho
    keys <- paste(selected, sprintf("%.17g", rho), sep = ":")
    stratum_map[[r]] <- keys
    sc_df <- deconvolution$scaling[[r]]
    valid <- sc_df$status == "ok"
    act_r <- deconvolution$activity[[r]]
    div_r <- sc_df$scale
    div_r[!is.finite(div_r) | div_r <= 0] <- 1
    noise_sd_r <- if (!is.null(sc_df$noise_sd)) sc_df$noise_sd else rep(NA_real_, ncol(act_r))
    need_nsd <- valid & (!is.finite(noise_sd_r) | noise_sd_r <= 0)
    if (any(need_nsd)) {
      noise_sd_r[need_nsd] <- div_r[need_nsd]
    }
    sig_eff_r <- if (isTRUE(restore_scale)) noise_sd_r else (noise_sd_r / div_r)
    sig_eff_r[!is.finite(sig_eff_r) | sig_eff_r <= 0] <- 1
    sig_eff_mat[r, ] <- sig_eff_r

    if (isTRUE(restore_scale)) {
      act_r[, valid] <- sweep(act_r[, valid, drop = FALSE], 2L, div_r[valid], "*")
    }

    dur_kernel <- if (readout_kernel == "linear" && source_model == "duration" && any(dur_target > 0)) {
      "boxcar"
    } else {
      readout_kernel
    }
    if (source_model == "duration" && any(duration > 0)) {
      S_all_0 <- build_epoch_readout_matrix(ev$onset + 0.5 * duration, 0, n, tr, offset,
        readout_kernel = dur_kernel, readout_width = pmax(w_src, duration))
    } else {
      S_all_0 <- build_epoch_readout_matrix(ev$onset, 0, n, tr, offset,
        readout_kernel = readout_kernel, readout_width = readout_width)
    }
    S_lags <- lapply(lags, function(lag_val) {
      build_epoch_readout_matrix(ev$onset[t_idx], lag_val, n, tr, offset,
        readout_kernel = readout_kernel, readout_width = readout_width)
    })
    S_lags_model <- if (source_model == "duration") {
      lapply(lags, function(lag_val) {
        build_epoch_readout_matrix(ev$onset[t_idx] + 0.5 * dur_target, lag_val, n, tr, offset,
          readout_kernel = dur_kernel, readout_width = pmax(w_src, dur_target))
      })
    } else {
      S_lags
    }

    for (li in seq_len(n_lags)) {
      raw_epochs[[li]][row_span, valid] <- S_lags[[li]] %*% act_r[, valid, drop = FALSE]
      demixed_epochs[[li]][row_span, valid] <- S_lags_model[[li]] %*% act_r[, valid, drop = FALSE]
    }

    for (key in unique(keys[valid])) {
      voxels <- which(keys == key & valid)
      k <- selected[voxels[1L]]
      ar <- rho[voxels[1L]]
      hrf_k <- fit$kernels[, k, drop = FALSE]
      operator <- hrf_inverse_operator(n, hrf_k[, 1L], rows, ar, z, settings)

      # Forward design for all events at their own onsets (for non-target leakage removal)
      dur_all_model <- duration
      if (source_model == "impulse") dur_all_model[t_idx] <- 0
      G_all_0 <- hrf_event_design_cpp(n, tr, offset, ev$onset, dur_all_model,
        seq_len(m_all), m_all, hrf_k)[[1L]][rows, , drop = FALSE]
      if (source_model == "duration" && duration_scale == "integral") {
        pos_t <- t_idx[dur_target > 0]
        if (length(pos_t)) {
          G_all_0[, pos_t] <- sweep(G_all_0[, pos_t, drop = FALSE], 2L,
            dur_target[dur_target > 0] / tr, "/")
        }
      }
      if (settings$scaling != "none") G_all_0 <- sweep(G_all_0, 2L, colMeans(G_all_0), "-")
      DG_all_0 <- operator %*% G_all_0

      L_other <- if (demix_other_events && length(o_idx) > 0L) {
        build_other_demix_operator(S_all_0, DG_all_0, ev, o_idx, t_idx, demix_ridge)
      } else {
        NULL
      }
      a_other_hat <- if (!is.null(L_other)) {
        L_other %*% act_r[, voxels, drop = FALSE]
      } else {
        NULL
      }

      # For duration source_model, estimate baseline firing rate alpha0_hat at lag = 0
      # so differential time-on-task baseline leakage is removed across all lags
      alpha0_hat <- NULL
      b_alpha0 <- NULL
      DG_target_0_sum <- NULL
      if (source_model == "duration") {
        S_0 <- build_epoch_readout_matrix(ev$onset[t_idx] + 0.5 * dur_target, 0, n, tr, offset,
          readout_kernel = dur_kernel, readout_width = pmax(w_src, dur_target))
        ok_0 <- rowSums(!is.finite(S_0)) == 0L
        G_0_unc <- hrf_event_design_cpp(n, tr, offset, ev$onset[t_idx], dur_target,
          seq_len(m_target), m_target, hrf_k)[[1L]][rows, , drop = FALSE]
        if (duration_scale == "integral" && any(dur_target > 0)) {
          d_pos0 <- dur_target > 0
          G_0_unc[, d_pos0] <- sweep(G_0_unc[, d_pos0, drop = FALSE], 2L, dur_target[d_pos0] / tr, "/")
        }
        G_0_cent <- if (settings$scaling != "none") {
          sweep(G_0_unc, 2L, colMeans(G_0_unc), "-")
        } else {
          G_0_unc
        }
        if (any(ok_0)) {
          DG_target_0_sum <- operator %*% rowSums(G_0_cent[, ok_0, drop = FALSE])
          Y_0 <- S_0[ok_0, , drop = FALSE] %*% act_r[, voxels, drop = FALSE]
          S_clean_0 <- S_0[ok_0, , drop = FALSE]
          if (!is.null(a_other_hat)) {
            R_other_0 <- S_0[ok_0, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
            Y_0 <- Y_0 - R_other_0 %*% a_other_hat
            S_clean_0 <- S_clean_0 - R_other_0 %*% L_other
          }
          demix_0 <- solve_duration_demix_operator(
            S_0[ok_0, , drop = FALSE], operator,
            G_0_unc[, ok_0, drop = FALSE], G_0_cent[, ok_0, drop = FALSE],
            settings$scaling, demix_ridge, gain_floor
          )
          alpha0_hat <- colMeans(demix_0$W %*% Y_0)
          b_alpha0 <- colMeans(demix_0$W %*% (S_clean_0 %*% operator))
        }
      }

      for (li in seq_len(n_lags)) {
        # Raw readout operator error variance
        S_raw_l <- S_lags[[li]]
        ok_raw <- rowSums(!is.finite(S_raw_l)) == 0L
        if (any(ok_raw)) {
          m_raw <- sum(ok_raw)
          B_raw <- S_raw_l[ok_raw, , drop = FALSE] %*% operator
          V_raw <- compute_ar1_operator_cov(B_raw, rows, ar, center_obs = center_obs)
          if (m_raw > 1L) {
            Q_raw <- diag(m_raw) - matrix(1 / m_raw, m_raw, m_raw)
            V_raw_cent <- (Q_raw %*% V_raw %*% Q_raw) / (1 - 1 / m_raw)
          } else {
            V_raw_cent <- V_raw
          }
          d_raw_unit <- sqrt(pmax(0, diag(V_raw_cent)))
          se_raw_epochs[[li]][row_span[ok_raw], voxels] <- outer(d_raw_unit, sig_eff_r[voxels])
        }

        S_tau <- S_lags_model[[li]]
        ok_tau <- rowSums(!is.finite(S_tau)) == 0L
        if (!any(ok_tau)) next
        m_ok <- sum(ok_tau)
        Y_clean <- demixed_epochs[[li]][row_span, voxels, drop = FALSE]
        S_clean_tau <- S_tau[ok_tau, , drop = FALSE]

        if (!is.null(a_other_hat)) {
          R_other_tau <- S_tau[ok_tau, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
          Y_clean[ok_tau, ] <- Y_clean[ok_tau, , drop = FALSE] - R_other_tau %*% a_other_hat
          S_clean_tau <- S_clean_tau - R_other_tau %*% L_other
        }
        B_clean_tau <- S_clean_tau %*% operator

        if (source_model == "duration") {
          shifted_onsets <- pmin(max(ev$onset), pmax(offset, ev$onset[t_idx] + lags[li]))
          G_tau_unc <- hrf_event_design_cpp(n, tr, offset, shifted_onsets, dur_target,
            seq_len(m_target), m_target, hrf_k)[[1L]][rows, , drop = FALSE]
          if (duration_scale == "integral" && any(dur_target > 0)) {
            d_pos <- dur_target > 0
            G_tau_unc[, d_pos] <- sweep(G_tau_unc[, d_pos, drop = FALSE], 2L, dur_target[d_pos] / tr, "/")
          }
          G_tau <- if (settings$scaling != "none") {
            sweep(G_tau_unc, 2L, colMeans(G_tau_unc), "-")
          } else {
            G_tau_unc
          }
          demix_tau <- solve_duration_demix_operator(
            S_tau[ok_tau, , drop = FALSE], operator,
            G_tau_unc[, ok_tau, drop = FALSE], G_tau[, ok_tau, drop = FALSE],
            settings$scaling, demix_ridge, gain_floor
          )
          K_norm <- demix_tau$K_norm
          g_vec <- diag(K_norm)
          mean_g <- mean(g_vec[is.finite(g_vec)])
          if (is.finite(mean_g)) {
            gain_accum[li, voxels] <- gain_accum[li, voxels] + mean_g
            gain_count[li, voxels] <- gain_count[li, voxels] + 1L
          }
          if (demix_target_trials) {
            Y_sub <- Y_clean[ok_tau, , drop = FALSE]
            B_sub_tau <- B_clean_tau
            if (!is.null(alpha0_hat) && !is.null(DG_target_0_sum) && abs(lags[li]) > 1e-8) {
              k_base_0 <- as.numeric(S_tau[ok_tau, , drop = FALSE] %*% DG_target_0_sum)
              k_base_tau <- rowSums(K_norm)
              delta_base <- k_base_0 - k_base_tau
              Y_sub <- Y_sub - outer(delta_base, alpha0_hat)
              if (!is.null(b_alpha0)) {
                B_sub_tau <- B_sub_tau - tcrossprod(delta_base, b_alpha0)
              }
            }
            A_hat <- demix_tau$W %*% Y_sub
            A_mean <- colMeans(A_hat)
            Y_clean[ok_tau, ] <- sweep(
              sweep(sweep(A_hat, 2L, A_mean, "-"), 1L, demix_tau$eff_gain_vec, "/"),
              2L, A_mean, "+"
            )
            Q_tau <- diag(m_ok) - matrix(1 / m_ok, m_ok, m_ok)
            G_gain_tau <- sweep(Q_tau, 1L, demix_tau$eff_gain_vec, "/") + (diag(m_ok) - Q_tau)
            B_demix_tau <- G_gain_tau %*% (demix_tau$W %*% B_sub_tau)
          } else {
            B_demix_tau <- B_clean_tau
          }
        } else {
          # Local target resolution matrix K_tau for unit sources at target onset + lag
          dur_tau <- if (source_model == "window") rep(w_src, m_target) else rep(0, m_target)
          shift_center <- if (source_model == "window") lags[li] - 0.5 * w_src else lags[li]
          shifted_onsets <- pmin(max(ev$onset), pmax(offset, ev$onset[t_idx] + shift_center))
          G_tau <- hrf_event_design_cpp(n, tr, offset, shifted_onsets, dur_tau,
            seq_len(m_target), m_target, hrf_k)[[1L]][rows, , drop = FALSE]
          if (source_model == "window" && w_src > 0) {
            G_tau <- G_tau / (w_src / tr)
          }
          if (settings$scaling != "none") G_tau <- sweep(G_tau, 2L, colMeans(G_tau), "-")
          K_norm <- S_tau[ok_tau, , drop = FALSE] %*% (operator %*% G_tau[, ok_tau, drop = FALSE])
          g_vec <- diag(K_norm)
          mean_g <- mean(g_vec[is.finite(g_vec)])
          if (is.finite(mean_g)) {
            gain_accum[li, voxels] <- gain_accum[li, voxels] + mean_g
            gain_count[li, voxels] <- gain_count[li, voxels] + 1L
          }

          if (demix_target_trials) {
            g_safe <- sign(ifelse(g_vec == 0, 1, g_vec)) * pmax(abs(g_vec), gain_floor)
            # Unit-diagonal target cross-talk matrix C_tau
            C_tau <- sweep(K_norm, 2L, g_safe, "/")
            lhs_c <- crossprod(C_tau) + demix_ridge * diag(nrow(C_tau))
            W_imp_tau <- sweep(solve(lhs_c, t(C_tau)), 1L, g_safe, "/")
            Y_clean[ok_tau, ] <- W_imp_tau %*% Y_clean[ok_tau, , drop = FALSE]
            B_demix_tau <- W_imp_tau %*% B_clean_tau
          } else {
            B_demix_tau <- B_clean_tau
          }
        }
        demixed_epochs[[li]][row_span, voxels] <- Y_clean

        # Compute dimensionless trial error covariance for stratum (r, key, li)
        V_demix <- compute_ar1_operator_cov(B_demix_tau, rows, ar, center_obs = center_obs)
        if (m_ok > 1L) {
          Q_tau <- diag(m_ok) - matrix(1 / m_ok, m_ok, m_ok)
          V_demix_cent <- (Q_tau %*% V_demix %*% Q_tau) / (1 - 1 / m_ok)
        } else {
          V_demix_cent <- V_demix
        }
        d_demix_unit <- sqrt(pmax(0, diag(V_demix_cent)))
        se_demixed_epochs[[li]][row_span[ok_tau], voxels] <- outer(d_demix_unit, sig_eff_r[voxels])
        if (m_ok > 1L) {
          num_ac1 <- V_demix_cent[cbind(2:m_ok, 1:(m_ok - 1L))]
          den_ac1 <- pmax(d_demix_unit[2:m_ok] * d_demix_unit[1:(m_ok - 1L)], .Machine$double.eps)
          ac1_vec <- c(0, pmin(0.999, pmax(-0.999, num_ac1 / den_ac1)))
        } else {
          ac1_vec <- 0
        }
        error_ac1_epochs[[li]][row_span[ok_tau], voxels] <- matrix(ac1_vec, nrow = m_ok, ncol = length(voxels))
        V_demix_blocks[[li]][[r]][[key]] <- list(ok_tau = ok_tau, V_cent = V_demix_cent)
      }
    }
  }

  has_gain <- gain_count > 0L
  gain_by_lag_voxel[has_gain] <- gain_accum[has_gain] / gain_count[has_gain]

  # Center covariates within run for pooled regression
  X_pool <- do.call(rbind, x_centered_list)
  coef_rows <- vector("list", n_lags * n_voxels * n_terms)
  row_k <- 0L

  for (li in seq_len(n_lags)) {
    Y_raw_lag <- raw_epochs[[li]]
    Y_demix_lag <- demixed_epochs[[li]]
    # Center epoch readouts within each run to eliminate run baseline offsets
    for (r in seq_len(n_runs)) {
      rs <- starts[r]:ends[r]
      Y_raw_lag[rs, ] <- sweep(Y_raw_lag[rs, , drop = FALSE], 2L,
        colMeans(Y_raw_lag[rs, , drop = FALSE], na.rm = TRUE), "-")
      Y_demix_lag[rs, ] <- sweep(Y_demix_lag[rs, , drop = FALSE], 2L,
        colMeans(Y_demix_lag[rs, , drop = FALSE], na.rm = TRUE), "-")
    }

    for (vi in seq_len(n_voxels)) {
      y_raw_v <- Y_raw_lag[, vi]
      y_dem_v <- Y_demix_lag[, vi]
      s_raw_v <- se_raw_epochs[[li]][, vi]
      s_dem_v <- se_demixed_epochs[[li]][, vi]
      s2_raw_v <- s_raw_v * s_raw_v
      s2_dem_v <- s_dem_v * s_dem_v

      ok_s_dem <- is.finite(s2_dem_v) & (s2_dem_v > 0)
      if (any(ok_s_dem)) {
        floor_op <- 0.10 * stats::median(s2_dem_v[ok_s_dem])
        w_op_vec <- 1 / (floor_op + s2_dem_v[ok_s_dem])
        pw_op_epochs[[li]][ok_s_dem, vi] <- w_op_vec / mean(w_op_vec)
      }
      ok_s_raw <- is.finite(s2_raw_v) & (s2_raw_v > 0)
      if (any(ok_s_raw)) {
        floor_raw <- 0.10 * stats::median(s2_raw_v[ok_s_raw])
        w_raw_vec <- 1 / (floor_raw + s2_raw_v[ok_s_raw])
        pw_raw_epochs[[li]][ok_s_raw, vi] <- w_raw_vec / mean(w_raw_vec)
      }

      ok_v <- is.finite(y_raw_v) & is.finite(y_dem_v) & (rowSums(!is.finite(X_pool)) == 0L)
      n_ok <- sum(ok_v)
      mean_se_v <- if (any(ok_s_dem)) mean(s_dem_v[ok_s_dem]) else NA_real_
      ac1_v <- error_ac1_epochs[[li]][, vi]
      mean_ac1_v <- if (any(is.finite(ac1_v))) mean(ac1_v[is.finite(ac1_v) & ac1_v != 0]) else NA_real_
      if (!is.finite(mean_ac1_v)) mean_ac1_v <- 0

      if (n_ok <= n_terms + 2L) {
        if (any(ok_s_dem)) {
          pw_epochs[[li]][ok_s_dem, vi] <- pw_op_epochs[[li]][ok_s_dem, vi]
        }
        for (ti in seq_len(n_terms)) {
          row_k <- row_k + 1L
          coef_rows[[row_k]] <- data.frame(
            voxel = fit$voxel_ids[vi], lag = lags[li], term = term_names[ti],
            beta_raw = NA_real_, beta_corrected = NA_real_,
            beta_unweighted = NA_real_, beta_weighted = NA_real_, beta_gls = NA_real_,
            se_corrected = NA_real_, se_unweighted = NA_real_, se_weighted = NA_real_, se_gls = NA_real_,
            r_raw = NA_real_, r_demixed = NA_real_,
            t_stat = NA_real_, p_value = NA_real_,
            t_unweighted = NA_real_, t_weighted = NA_real_, t_gls = NA_real_,
            df = n_ok - n_terms,
            mean_gain = gain_by_lag_voxel[li, vi],
            mean_se_trial = mean_se_v,
            mean_error_ac1 = mean_ac1_v,
            stringsAsFactors = FALSE
          )
        }
        next
      }

      X_ok <- X_pool[ok_v, , drop = FALSE]
      xtx_raw <- crossprod(X_ok) + 1e-8 * diag(n_terms)
      xtx_inv_raw <- solve(xtx_raw)
      beta_raw_v <- as.numeric(xtx_inv_raw %*% crossprod(X_ok, y_raw_v[ok_v]))
      beta_unw_v <- as.numeric(xtx_inv_raw %*% crossprod(X_ok, y_dem_v[ok_v]))

      resid_unw <- y_dem_v[ok_v] - as.numeric(X_ok %*% beta_unw_v)
      df_v <- max(1L, n_ok - n_terms - n_runs)
      sigma2_unw <- sum(resid_unw^2) / df_v
      se_unw_naive <- sqrt(pmax(diag(xtx_inv_raw) * sigma2_unw, .Machine$double.eps))

      # Assemble per-run operator covariance blocks and estimate biological variance tau2_bio
      Omega_r_list <- vector("list", n_runs)
      rs_ok_list <- vector("list", n_runs)
      xOx_sum <- matrix(0, n_terms, n_terms)
      tr_Omega <- 0
      for (r in seq_len(n_runs)) {
        rs <- starts[r]:ends[r]
        ok_r_local <- which(ok_v[rs])
        if (!length(ok_r_local)) next
        rs_ok_list[[r]] <- rs[ok_r_local]
        key_rv <- if (!is.null(stratum_map[[r]])) stratum_map[[r]][vi] else NULL
        blk_r <- if (!is.null(key_rv)) V_demix_blocks[[li]][[r]][[key_rv]] else NULL
        if (!is.null(blk_r)) {
          ok_tau_idx <- which(blk_r$ok_tau)
          pos_in_blk <- match(ok_r_local, ok_tau_idx)
          valid_pos <- !is.na(pos_in_blk)
          if (all(valid_pos)) {
            Omega_r <- (sig_eff_mat[r, vi]^2) * blk_r$V_cent[pos_in_blk, pos_in_blk, drop = FALSE]
          } else {
            s2_sub <- s2_dem_v[rs[ok_r_local]]
            s2_sub[!is.finite(s2_sub)] <- 1
            Omega_r <- diag(s2_sub, length(ok_r_local))
          }
        } else {
          s2_sub <- s2_dem_v[rs[ok_r_local]]
          s2_sub[!is.finite(s2_sub)] <- 1
          Omega_r <- diag(s2_sub, length(ok_r_local))
        }
        Omega_r_list[[r]] <- Omega_r
        tr_Omega <- tr_Omega + sum(diag(Omega_r))
        Xr_ok <- X_pool[rs[ok_r_local], , drop = FALSE]
        xOx_sum <- xOx_sum + crossprod(Xr_ok, Omega_r %*% Xr_ok)
      }
      tr_proj_Omega <- max(0, tr_Omega - sum(diag(xtx_inv_raw %*% xOx_sum)))
      tau2_bio <- max(0, (sum(resid_unw^2) - tr_proj_Omega) / df_v)
      s2_ok <- s2_dem_v[ok_v]
      s2_med <- if (any(is.finite(s2_ok) & s2_ok > 0)) stats::median(s2_ok[is.finite(s2_ok) & s2_ok > 0]) else 1
      tau2_eff <- max(tau2_bio, 0.10 * s2_med)

      if (any(ok_s_dem)) {
        w_eb_all <- 1 / (tau2_eff + s2_dem_v[ok_s_dem])
        pw_epochs[[li]][ok_s_dem, vi] <- w_eb_all / mean(w_eb_all)
      } else {
        pw_epochs[[li]][ok_v, vi] <- 1
      }
      w_v_all <- pw_epochs[[li]][, vi]
      w_v_all[!is.finite(w_v_all) | w_v_all <= 0] <- 1

      # Precision-weighted WLS and within-run ZCA-whitened GLS
      X_w_ok <- matrix(0, n_ok, n_terms)
      y_w_ok <- numeric(n_ok)
      X_gls_ok <- matrix(0, n_ok, n_terms)
      y_gls_ok <- numeric(n_ok)
      w_ok <- numeric(n_ok)
      Sigma_tot_list <- vector("list", n_runs)
      pos_cursor <- 0L

      for (r in seq_len(n_runs)) {
        rs_ok <- rs_ok_list[[r]]
        m_r <- length(rs_ok)
        if (!m_r) next
        idx_span <- (pos_cursor + 1L):(pos_cursor + m_r)
        pos_cursor <- pos_cursor + m_r

        wr <- w_v_all[rs_ok]
        w_ok[idx_span] <- wr

        yr_raw <- demixed_epochs[[li]][rs_ok, vi]
        Xr_raw <- X_raw_pool[rs_ok, , drop = FALSE]
        mu_y_w <- sum(wr * yr_raw) / sum(wr)
        mu_X_w <- colSums(Xr_raw * wr) / sum(wr)
        yr_w <- yr_raw - mu_y_w
        Xr_w <- sweep(Xr_raw, 2L, mu_X_w, "-")
        y_w_ok[idx_span] <- yr_w
        X_w_ok[idx_span, ] <- Xr_w

        Omega_r <- Omega_r_list[[r]]
        tau2_zca <- max(tau2_eff, 0.35 * sigma2_unw)
        Sigma_tot_r <- Omega_r + diag(tau2_zca, m_r)
        Sigma_tot_list[[r]] <- list(idx = idx_span, Sigma = Sigma_tot_r, wr = wr)

        Phi_r <- compute_zca_whitening_matrix(Sigma_tot_r, eig_floor = 0.05)
        yr_cent <- yr_raw - mean(yr_raw)
        Xr_cent <- sweep(Xr_raw, 2L, colMeans(Xr_raw), "-")
        yr_zca <- as.numeric(Phi_r %*% yr_cent)
        Xr_zca <- Phi_r %*% Xr_cent
        if (isTRUE(whiten_panel)) {
          whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + yr_zca
          for (ti in seq_len(n_terms)) {
            whitened_x_epochs[[li]][[ti]][rs_ok, vi] <- mean(Xr_raw[, ti]) + Xr_zca[, ti]
          }
        }
        mu_y_gls <- sum(wr * yr_zca) / sum(wr)
        mu_X_gls <- colSums(Xr_zca * wr) / sum(wr)
        y_gls_ok[idx_span] <- yr_zca - mu_y_gls
        X_gls_ok[idx_span, ] <- sweep(Xr_zca, 2L, mu_X_gls, "-")
      }

      w_ok <- w_ok / mean(w_ok)
      WX_w <- X_w_ok * w_ok
      xwx <- crossprod(X_w_ok, WX_w) + 1e-8 * diag(n_terms)
      xwx_inv <- solve(xwx)
      beta_pw_v <- as.numeric(xwx_inv %*% crossprod(WX_w, y_w_ok))
      resid_pw <- y_w_ok - as.numeric(X_w_ok %*% beta_pw_v)

      # Covariance-propagated sandwich SE for precision-weighted estimator
      meat_pw <- matrix(0, n_terms, n_terms)
      tr_W_Sigma <- 0
      for (r in seq_len(n_runs)) {
        st_r <- Sigma_tot_list[[r]]
        if (is.null(st_r)) next
        idx_r <- st_r$idx
        WX_r <- WX_w[idx_r, , drop = FALSE]
        Sigma_r <- st_r$Sigma
        wr <- w_ok[idx_r]
        meat_r <- crossprod(WX_r, Sigma_r %*% WX_r)
        meat_pw <- meat_pw + meat_r
        tr_W_Sigma <- tr_W_Sigma + sum(wr * diag(Sigma_r)) - sum(diag(xwx_inv %*% meat_r))
      }
      scale_pw <- max(0.25, sum(w_ok * (resid_pw^2)) / max(tr_W_Sigma, 1e-8))
      cov_beta_pw <- scale_pw * (xwx_inv %*% meat_pw %*% xwx_inv)
      se_pw_v <- sqrt(pmax(diag(cov_beta_pw), .Machine$double.eps))

      # ZCA-whitened block-GLS estimator
      WX_gls <- X_gls_ok * w_ok
      xwx_gls <- crossprod(X_gls_ok, WX_gls) + 1e-8 * diag(n_terms)
      xwx_gls_inv <- solve(xwx_gls)
      beta_gls_v <- as.numeric(xwx_gls_inv %*% crossprod(WX_gls, y_gls_ok))
      resid_gls <- y_gls_ok - as.numeric(X_gls_ok %*% beta_gls_v)
      sigma2_gls <- sum(w_ok * (resid_gls^2)) / df_v
      se_gls_v <- sqrt(pmax(diag(xwx_gls_inv) * sigma2_gls, .Machine$double.eps))

      t_unw_vec <- beta_unw_v / se_unw_naive
      t_pw_vec <- beta_pw_v / se_pw_v
      t_gls_vec <- beta_gls_v / se_gls_v

      if (identical(weighting, "precision")) {
        beta_corr_v <- beta_pw_v
        se_corr_v <- se_pw_v
        t_corr_vec <- t_pw_vec
      } else if (identical(weighting, "gls")) {
        beta_corr_v <- beta_gls_v
        se_corr_v <- se_gls_v
        t_corr_vec <- t_gls_vec
      } else {
        beta_corr_v <- beta_unw_v
        se_corr_v <- se_unw_naive
        t_corr_vec <- t_unw_vec
      }
      p_corr_vec <- 2 * stats::pt(abs(t_corr_vec), df = df_v, lower.tail = FALSE)

      for (ti in seq_len(n_terms)) {
        xk <- X_ok[, ti]
        r_raw_val <- if (stats::sd(xk) > 0 && stats::sd(y_raw_v[ok_v]) > 0) {
          stats::cor(xk, y_raw_v[ok_v])
        } else NA_real_
        r_dem_val <- if (stats::sd(xk) > 0 && stats::sd(y_dem_v[ok_v]) > 0) {
          stats::cor(xk, y_dem_v[ok_v])
        } else NA_real_
        row_k <- row_k + 1L
        coef_rows[[row_k]] <- data.frame(
          voxel = fit$voxel_ids[vi], lag = lags[li], term = term_names[ti],
          beta_raw = beta_raw_v[ti], beta_corrected = beta_corr_v[ti],
          beta_unweighted = beta_unw_v[ti], beta_weighted = beta_pw_v[ti], beta_gls = beta_gls_v[ti],
          se_corrected = se_corr_v[ti], se_unweighted = se_unw_naive[ti],
          se_weighted = se_pw_v[ti], se_gls = se_gls_v[ti],
          r_raw = r_raw_val, r_demixed = r_dem_val,
          t_stat = t_corr_vec[ti], p_value = p_corr_vec[ti],
          t_unweighted = t_unw_vec[ti], t_weighted = t_pw_vec[ti], t_gls = t_gls_vec[ti],
          df = df_v,
          mean_gain = gain_by_lag_voxel[li, vi],
          mean_se_trial = mean_se_v,
          mean_error_ac1 = mean_ac1_v,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  coef_df <- do.call(rbind, coef_rows)
  fmean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  summary_rows <- list()
  for (li in seq_len(n_lags)) {
    for (ti in seq_len(n_terms)) {
      sub_df <- coef_df[coef_df$lag == lags[li] & coef_df$term == term_names[ti], , drop = FALSE]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        lag = lags[li], term = term_names[ti],
        mean_beta_raw = fmean(sub_df$beta_raw),
        mean_beta_corrected = fmean(sub_df$beta_corrected),
        mean_beta_unweighted = fmean(sub_df$beta_unweighted),
        mean_beta_weighted = fmean(sub_df$beta_weighted),
        mean_beta_gls = fmean(sub_df$beta_gls),
        mean_r_raw = fmean(sub_df$r_raw),
        mean_r_demixed = fmean(sub_df$r_demixed),
        mean_gain = fmean(sub_df$mean_gain),
        mean_se_trial = fmean(sub_df$mean_se_trial),
        mean_error_ac1 = fmean(sub_df$mean_error_ac1),
        n_voxels = sum(is.finite(sub_df$beta_corrected)),
        stringsAsFactors = FALSE
      )
    }
  }

  panel_dt <- NULL
  if (isTRUE(return_panel)) {
    run_vec <- rep(fit$run_ids, times = target_counts)
    trial_vec <- unlist(lapply(target_counts, seq_len), use.names = FALSE)
    evt_idx_vec <- unlist(target_idx_list, use.names = FALSE)
    onset_vec <- unlist(lapply(seq_len(n_runs), function(r) events[[r]]$onset[target_idx_list[[r]]]), use.names = FALSE)
    cov_full <- as.data.frame(do.call(rbind, cov_list))

    panel_list <- vector("list", n_lags)
    for (li in seq_len(n_lags)) {
      base_meta <- data.table::data.table(
        run = rep(run_vec, times = n_voxels),
        trial = rep(trial_vec, times = n_voxels),
        event_index = rep(evt_idx_vec, times = n_voxels),
        onset = rep(onset_vec, times = n_voxels),
        evt_time = lags[li],
        voxel = rep(fit$voxel_ids, each = total_trials),
        decon_raw = as.numeric(raw_epochs[[li]]),
        decon_demixed = as.numeric(demixed_epochs[[li]]),
        se_raw = as.numeric(se_raw_epochs[[li]]),
        se_demixed = as.numeric(se_demixed_epochs[[li]]),
        precision_weight = as.numeric(pw_epochs[[li]]),
        precision_weight_op = as.numeric(pw_op_epochs[[li]]),
        precision_weight_raw = as.numeric(pw_raw_epochs[[li]]),
        error_ac1 = as.numeric(error_ac1_epochs[[li]])
      )
      for (cn in colnames(cov_full)) {
        base_meta[[cn]] <- rep(cov_full[[cn]], times = n_voxels)
      }
      if (isTRUE(whiten_panel)) {
        base_meta[["decon_demixed_whitened"]] <- as.numeric(whitened_y_epochs[[li]])
        for (ti in seq_len(n_terms)) {
          cn_w <- paste0(term_names[ti], "_whitened")
          base_meta[[cn_w]] <- as.numeric(whitened_x_epochs[[li]][[ti]])
        }
      }
      panel_list[[li]] <- base_meta
    }
    panel_dt <- data.table::rbindlist(panel_list)
  }

  list(
    coefficients = coef_df,
    summary = do.call(rbind, summary_rows),
    epochs = if (return_epochs) {
      list(
        raw = raw_epochs,
        demixed = demixed_epochs,
        se_raw = se_raw_epochs,
        se_demixed = se_demixed_epochs,
        precision_weight = pw_epochs,
        precision_weight_op = pw_op_epochs,
        error_ac1 = error_ac1_epochs
      )
    } else NULL,
    panel = panel_dt,
    settings = list(
      target_condition = target_condition,
      lags = lags,
      readout_kernel = readout_kernel,
      readout_width = readout_width,
      demix_other_events = demix_other_events,
      demix_target_trials = demix_target_trials,
      demix_ridge = demix_ridge,
      gain_floor = gain_floor,
      restore_scale = restore_scale,
      source_model = source_model,
      duration_scale = duration_scale,
      weighting = weighting,
      whiten_panel = whiten_panel
    )
  )
}
