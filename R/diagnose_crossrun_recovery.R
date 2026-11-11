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
    # Identify chronologically preceding and succeeding events (across all condition types)
    # so we can measure immediate neighbor cross-talk in the resolution matrix.
    previous <- ifelse(position > 1L, order[pmax(1L, position - 1L)], NA_integer_)
    next_event <- ifelse(position < m, order[pmin(m, position + 1L)], NA_integer_)
    # Build linear readout matrix S (m_events x n_time) using two-bin fractional
    # bracket interpolation at each event onset + lag.
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
    # Construct the unblurred neural design B (using a delta kernel [1, 0]) and
    # compute ideal = S * B so that natural overlap of neural readout windows is
    # separated from HRF + regularization blur (excess = S * D * H * B - S * B).
    neural <- hrf_event_design_cpp(n, tr, offset, ev$onset, duration,
      seq_len(m), m, matrix(c(1, 0), ncol = 1L))[[1L]]
    ideal <- readout %*% neural
    ideal_diagonal <- diag(ideal)
    mask <- if (settings$observation_mask == "fit") fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    if (length(rows) < 3L) next
    z <- if (settings$nuisance_mode == "joint") cbind(1, seq(-1, 1, length.out = n),
      fit$specification$nuisance[[r]]) else matrix(numeric(), n, 0L)
    strat <- hrf_stratum_kernels(fit, r)
    selected <- fit$selections[[r]]$kernel_index
    rho <- strat$rho
    # Group voxels by their (selected_hrf_kernel / weight_vector, ar1_rho) stratum so the
    # inverse operator D is constructed and multiplied only once per unique stratum.
    keys <- strat$keys
    valid <- deconvolution$scaling[[r]]$status == "ok"
    mappings[[r]] <- data.frame(run = fit$run_ids[r], voxel = fit$voxel_ids,
      stratum = keys, usable = valid)
    for (key in unique(keys[valid])) {
      voxels <- which(keys == key & valid); k <- selected[voxels[1L]]; ar <- rho[voxels[1L]]
      ar_op <- if (!is.null(strat$get_ar)) strat$get_ar(voxels[1L]) else ar
      hrf_k <- strat$get_kernel(voxels[1L])
      operator <- hrf_inverse_operator(n, as.numeric(hrf_k), rows, ar_op, z, settings)
      readout_operator <- readout %*% operator
      source <- if (is.null(sources)) hrf_k else sources
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

# Build a trial-by-time epoch readout matrix S(tau) of dimension (m_events x n_time).
# For each event onset t_i and relative lag tau, S(tau) maps the continuous
# deconvolved time series (length n_time, sampled every TR seconds starting at
# time_offset) onto a single trial-level readout value at t_i + tau.
# - "linear" performs two-point linear bracket interpolation between adjacent TRs.
# - "boxcar", "triangle", and "gaussian" evaluate a 9-point symmetric quadrature
#   over [t_i + tau - 0.5*width, t_i + tau + 0.5*width], attenuating high-frequency
#   deconvolution noise while preserving local event timing.
# - When readout_width is a vector of length(onsets) (used by source_model = "duration"),
#   each trial i receives its own width_scale[i] = max(w_src, duration[i]).
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
    # Fractional grid index for quadrature node q of each trial
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
  # Normalize row weights for trials partially supported near run boundaries;
  # trials completely outside the run grid receive NA.
  valid <- weight_sum > 0
  if (any(valid)) {
    readout[valid, ] <- sweep(readout[valid, , drop = FALSE], 1L, weight_sum[valid], "/")
  }
  readout[!valid, ] <- NA_real_
  readout
}

# Construct the target-trial cross-talk de-mixing and gain-correction operator W(tau)
# for a given lag tau.
# Given readout matrix S_tau (m_ok x n_time), inverse deconvolution operator D (n_time x n_obs),
# uncentered forward trial design G_tau_unc (n_obs x m_ok), and column-centered
# design G_tau_cent = G_tau_unc - 1 * mu_G^T:
# 1. K_unc = S_tau * D * G_tau_unc is the uncentered trial-by-trial resolution matrix,
#    whose diagonal g_unc[i] is trial i's isolated amplitude gain and whose off-diagonals
#    measure adjacent-trial overlap leakage.
# 2. Row-normalizing R_unc = diag(g_safe)^{-1} K_unc places all trials on unit diagonal
#    even when trial durations vary (where g_unc[i] scales with duration d_i).
# 3. Penalizing only centered trial fluctuations via Q = I - (1/m)*1*1^T leaves a
#    constant baseline firing rate alpha_0 * 1 unattenuated (since Q * 1 = 0),
#    preventing duration-dependent shrinkage of the task baseline from confounding
#    covariates correlated with reaction time (RT).
# 4. When input BOLD was centered (scaling != "none"), G_tau_cent differs from
#    G_tau_unc by the rank-1 matrix -1 * mu_G^T. A Sherman-Morrison rank-1 update
#    adjusts W_unc into W_full so that constant baseline activity still maps to a
#    constant across trials of unequal duration.
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
    # Sherman-Morrison rank-1 adjustment for column-demeaned forward design G_tau_cent
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
  # Effective post-demixing gain in the centered trial subspace Q
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

# Propagate observation-level AR(1) noise covariance Sigma_AR1(rho) through a
# linear trial readout/de-mixing operator B (M_trials x m_obs) in O(m_obs * M)
# operations without forming a dense m_obs x m_obs covariance matrix.
# Because an AR(1) process on retained observation indices `rows` has a Markovian
# forward innovation representation e_{rows[i]} = a_i * e_{rows[i-1]} + s_i * z_i
# with a_i = rho^(rows[i] - rows[i-1]) and s_i = sqrt(1 - a_i^2), the backward
# accumulation V[i, ] = B[, i]^T + a_{i+1} * V[i+1, ] satisfies B * Sigma * B^T = U^T * U
# where U[i, ] = s_i * V[i, ]. This yields the exact M x M trial error covariance.
compute_ar1_operator_factor <- function(B, rows, rho, center_obs = TRUE) {
  if (length(rows) > 0L && ncol(B) > length(rows) && max(rows) <= ncol(B)) {
    B <- B[, rows, drop = FALSE]
  }
  m <- ncol(B)
  M <- nrow(B)
  if (m == 0L || M == 0L) return(matrix(0, m, M))
  B_eff <- if (isTRUE(center_obs)) {
    sweep(B, 1L, rowMeans(B), "-")
  } else {
    B
  }
  rho1 <- if (length(rho)) rho[1L] else 0.0
  has_ar1 <- is.finite(rho1) && abs(rho1) >= 1e-10
  has_psi <- m >= 2L && length(rho) >= 2L && any(is.finite(rho[-1L]) & abs(rho[-1L]) >= 1e-10)
  if ((!has_ar1 && !has_psi) || m == 1L) {
    return(t(B_eff))
  }
  Bt <- t(B_eff)
  if (has_psi) {
    is_contig <- all(diff(rows) == 1L)
    lag_span <- if (is_contig) max(8L, length(rho) + 2L) else max(8L, max(rows) - min(rows))
    info <- hrf_ar_cascade_info(rho, max_lag = lag_span)
    phi <- info$phi
    p <- length(phi)
    if (p == 0L) return(Bt)
    if (is_contig && m > p) {
      Gamma_p <- matrix(info$acf[abs(outer(seq_len(p), seq_len(p), "-")) + 1L], p, p)
      L_p <- tryCatch(t(chol(Gamma_p)), error = function(e) NULL)
      if (!is.null(L_p) && is.finite(info$sig_u2) && info$sig_u2 > 1e-8) {
        sig_u <- sqrt(info$sig_u2)
        W_v <- matrix(0.0, m, M)
        for (i in m:1L) {
          acc <- Bt[i, ]
          for (k in seq_len(p)) {
            if (i + k <= m && i + k > p) {
              acc <- acc + phi[k] * W_v[i + k, ]
            }
          }
          W_v[i, ] <- acc
        }
        U <- matrix(0.0, m, M)
        U[(p + 1L):m, ] <- sig_u * W_v[(p + 1L):m, , drop = FALSE]
        U[seq_len(p), ] <- crossprod(L_p, W_v[seq_len(p), , drop = FALSE])
        return(U)
      }
    }
    idx <- abs(outer(rows, rows, "-")) + 1L
    Sigma_sub <- matrix(info$acf[idx], m, m)
    R_chol <- tryCatch(chol(Sigma_sub), error = function(e) {
      diag(Sigma_sub) <- diag(Sigma_sub) + 1e-6
      chol(Sigma_sub)
    })
    return(R_chol %*% Bt)
  }
  V <- matrix(0, m, M)
  V[m, ] <- Bt[m, ]
  a_vec <- c(0, rho1^(rows[-1L] - rows[-m]))
  for (i in (m - 1L):1L) {
    V[i, ] <- Bt[i, ] + a_vec[i + 1L] * V[i + 1L, ]
  }
  s_vec <- sqrt(pmax(0, 1 - a_vec * a_vec))
  sweep(V, 1L, s_vec, "*")
}

compute_ar1_operator_cov <- function(B, rows, rho, center_obs = TRUE) {
  m <- ncol(B)
  M <- nrow(B)
  if (m == 0L || M == 0L) return(matrix(0, M, M))
  U <- compute_ar1_operator_factor(B, rows, rho, center_obs = center_obs)
  crossprod(U)
}

# Construct linear operators that estimate non-target event amplitudes (L_other)
# and joint lag-0 target amplitudes (L_target_0) from deconvolved activity.
# In rapid multi-phase trials (e.g., Cue -> Decision -> Feedback), preceding or
# following non-target events leak into the target epoch through the HRF/inverse
# resolution matrix R_all_0 = S_all_0 * D * G_all_0.
# 1. First, condition-level mean amplitudes A_mean_op are estimated across all
#    conditions simultaneously so condition baselines are properly separated.
# 2. Next, residual trial-to-trial fluctuations Res_all_op are inverted using a
#    column-normalized, second-order Neumann-refined ridge inverse (W_2) to
#    capture trial-specific amplitude variation in non-target events (and, for
#    impulse models at lag 0, joint target-trial amplitudes).
build_other_demix_operator <- function(S_all_0, DG_all_0, ev, o_idx, t_idx, demix_ridge,
    source_model = "impulse", demix_target_trials = TRUE) {
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
  if (!length(ok_o) || !isTRUE(demix_target_trials) || demix_ridge >= 0.03) {
    return(list(
      L_other = A_mean_op[o_idx, , drop = FALSE],
      L_target_0 = NULL
    ))
  }
  Res_all_op <- S_ok - R_all_0 %*% A_mean_op
  diag_all <- pmax(diag(R_all_0), 0.05)
  C_all <- sweep(R_all_0, 2L, diag_all, "/")
  r_eff <- min(max(demix_ridge, 0.025), 0.035)
  W_1 <- solve(crossprod(C_all) + r_eff * diag(ncol(C_all)), t(C_all))
  # Second-order iterated Tikhonov / Neumann refinement reduces ridge shrinkage bias
  W_2 <- W_1 + 0.80 * W_1 %*% (diag(nrow(C_all)) - C_all %*% W_1)
  g_eff <- pmax(diag(W_2 %*% C_all), 0.10)
  Delta_all_op <- sweep(W_2 %*% Res_all_op, 1L, g_eff * diag_all, "/")
  alpha_delta <- if (identical(source_model, "impulse")) 1.0 else 0.0
  A_joint_op <- A_mean_op + alpha_delta * Delta_all_op
  list(
    L_other = A_joint_op[o_idx, , drop = FALSE],
    L_target_0 = if (identical(source_model, "impulse")) A_joint_op[t_idx, , drop = FALSE] else NULL
  )
}

# Compute a variance-preserving Zero-Phase Component Analysis (ZCA) whitening
# matrix Phi for a within-run trial covariance matrix Sigma_tot (M x M).
# Let D = diag(sqrt(diag(Sigma_tot))) and R_tot = D^{-1} Sigma_tot D^{-1} be the
# trial correlation matrix. Then Phi = D * R_tot^{-1/2} * D^{-1} (normalized to
# unit mean diagonal) removes off-diagonal trial-to-trial error correlations
# while preserving each trial's marginal error variance diag(Sigma_tot).
# This allows downstream mixed models (`mixed_by(..., weights = "precision_weight")`)
# to pair ZCA-decorrelated trial observations with diagonal precision weights.
compute_zca_whitening_matrix <- function(Sigma_tot, eig_floor = 0.05, offdiag_scale = 1.0) {
  M <- nrow(Sigma_tot)
  if (M <= 1L) return(diag(M))
  d_tot <- sqrt(pmax(diag(Sigma_tot), .Machine$double.eps))
  R_tot <- sweep(sweep(Sigma_tot, 1L, d_tot, "/"), 2L, d_tot, "/")
  R_tot <- 0.5 * (R_tot + t(R_tot))
  if (offdiag_scale < 1.0) {
    R_tot <- (1.0 - offdiag_scale) * diag(M) + offdiag_scale * R_tot
  }
  eig <- eigen(R_tot, symmetric = TRUE)
  cut_null <- max(eig_floor, 0.35)
  inv_sqrt_vals <- ifelse(eig$values < cut_null, sqrt(pmax(eig$values, 0)), 1 / sqrt(pmax(eig$values, eig_floor)))
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
#' relative to a target event while removing preceding/neighboring event leakage,
#' correcting regularization-induced amplitude attenuation, and propagating
#' per-trial operator precision weights and within-run error covariance.
#'
#' @param deconvolution Result of \code{deconvolve_crossrun_hrf()}.
#' @param covariates Optional numeric vector, matrix, or data.frame of
#'   trial-level predictors for target trials (concatenated across runs or
#'   supplied as a run-named list matching \code{deconvolution$hrf_fit$run_ids}).
#'   When \code{NULL} (default), performs covariate-agnostic epoch readout,
#'   de-mixing, biological variance estimation (\code{tau2_mode = "marginal"}),
#'   Empirical-Bayes precision weighting, and optional within-run ZCA whitening
#'   (\code{whiten_panel = TRUE}), returning the trial panel (\code{return_panel}
#'   defaults to \code{TRUE} when \code{covariates = NULL}) for downstream
#'   two-stage \code{\link{mixed_by}} behavioral modeling or inter-region
#'   connectivity (\code{R2 ~ R1_seed * x}).
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
#'   contributions at their own onsets and subtract their operator-predicted
#'   leakage into each target epoch lag. When \code{demix_ridge >= 0.03}
#'   (default \code{0.05}), uses a condition-mean non-target de-mixing operator
#'   that removes preceding cue/event bleed without inflating target-trial
#'   variance; when \code{demix_ridge < 0.03} (e.g. \code{0.02}), activates a
#'   2nd-order iterated Tikhonov joint \eqn{(T, O)} operator to de-mix strongly
#'   trial-varying non-target confounds. Default: \code{TRUE}.
#' @param demix_target_trials Logical; whether to use the stratum resolution
#'   matrix at each lag to de-mix adjacent target-trial overlap and restore
#'   regularization-attenuated amplitude gain via 2nd-order iterated Tikhonov
#'   inversion (\code{solve_duration_demix_operator()}). Default: \code{TRUE}.
#' @param demix_within_trial Logical; whether to apply an additional multi-lag
#'   point-spread de-mixing step across the requested \code{lags} (when at least
#'   3 lags are evaluated) to sharpen temporal localization across lags within
#'   each target trial. Recommended (\code{TRUE}) when estimating within-trial
#'   temporal waveforms \eqn{\beta(\tau)} across a multi-lag grid. Default:
#'   \code{FALSE}.
#' @param demix_ridge Positive ridge stabilization parameter for target and
#'   non-target cross-talk inversion. Values \code{>= 0.03} (default \code{0.05})
#'   use condition-mean non-target de-mixing; values \code{< 0.03} activate
#'   joint trial-varying non-target de-mixing. Default: \code{0.05}.
#' @param within_trial_ridge Positive ridge stabilization parameter for
#'   multi-lag within-trial point-spread inversion when
#'   \code{demix_within_trial = TRUE}. Default: \code{0.06}.
#' @param gain_floor Minimum absolute operator gain used when correcting
#'   amplitude slopes. Default: \code{0.05}.
#' @param restore_scale Logical; if \code{TRUE}, multiplies each run and voxel by
#'   its recorded \code{deconvolution$scaling[[r]]$scale} divisor so that
#'   \code{beta_raw} and \code{beta_corrected} (and panel readouts) are reported
#'   in unscaled physical input units rather than standardized
#'   (\code{zscore}/\code{noise_sd}) units. Strongly recommended (\code{TRUE})
#'   when pooling across multiple runs in \code{\link{mixed_by}}. Default:
#'   \code{FALSE}.
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
#'   errors and run noise scale), \code{"gls"} (null-mode-damped,
#'   Bartlett-adaptive within-run ZCA correlation whitening plus precision
#'   weighting accounting for off-diagonal trial-to-trial error correlations),
#'   or \code{"none"} (unweighted OLS). All three slope estimates
#'   (\code{beta_unweighted}, \code{beta_weighted}, \code{beta_gls}) and their
#'   standard errors are always included in \code{coefficients}.
#' @param tau2_mode Estimator for biological trial-to-trial variance
#'   \eqn{\hat{\tau}^2_{\mathrm{bio}}} used in Empirical-Bayes trial precision
#'   weights and within-run ZCA whitening: \code{"auto"} (default; uses
#'   \code{"marginal"} when \code{covariates = NULL} and an Empirical-Bayes
#'   blend of residual and marginal estimators when covariates are supplied),
#'   \code{"residual"} (method-of-moments estimator from OLS covariate regression
#'   residuals minus projected operator variance \eqn{\mathrm{tr}((\mathbf{I}-\mathbf{P}_X)\boldsymbol{\Omega})}),
#'   or \code{"marginal"} (covariate-agnostic estimator combining run-centered
#'   marginal variance, lag-1 trial difference variance, and robust MAD scale
#'   minus analytical operator variance expectations).
#' @param aggregate_by_parcel Logical; if \code{TRUE}, aggregates deconvolved
#'   signals across voxels within each spatial parcel (using inverse-variance
#'   voxel precision weights and empirical spatial noise correlation propagation)
#'   prior to Stage-2 epoch de-mixing and Stage-3 regression/ZCA whitening.
#'   Replaces voxel-by-voxel Stage-3/4 fitting with parcel-level inference while
#'   calibrating \code{se_demixed}, \code{tau2_bio}, \code{precision_weight},
#'   and ZCA whitening directly on the high-SNR parcel series. Default:
#'   \code{FALSE}.
#' @param parcel_groups Optional character or factor vector of parcel/region
#'   labels for each voxel (either unnamed of length \code{n_voxels} or named by
#'   \code{voxel_ids}). Defaults to
#'   \code{deconvolution$hrf_fit$specification$spatial_groups} when
#'   \code{aggregate_by_parcel = TRUE}. Supplying a non-\code{NULL}
#'   \code{parcel_groups} automatically enables \code{aggregate_by_parcel = TRUE}.
#' @param whiten_panel Logical; if \code{TRUE} and \code{return_panel = TRUE},
#'   adds within-run null-mode-damped, Bartlett-adaptive ZCA-whitened outcome
#'   (\code{decon_demixed_whitened}) and covariate (\code{<var>_whitened})
#'   columns to \code{panel} so downstream
#'   \code{mixed_by(..., weights = "precision_weight")} models account for both
#'   heteroscedastic trial precision and off-diagonal within-run trial error
#'   correlations. Recommended (\code{TRUE}) for group mixed-effects inference.
#'   Default: \code{FALSE}.
#' @param return_epochs Logical; whether to include raw and de-mixed trial-by-voxel
#'   (or trial-by-parcel) epoch matrices per lag. Default: \code{FALSE}.
#' @param return_panel Logical; whether to return a tidy MEDUSA-compatible
#'   \code{data.table} (\code{panel}) with columns \code{run}, \code{trial},
#'   \code{event_index}, \code{onset}, \code{evt_time}, \code{voxel} (and
#'   \code{parcel} when \code{aggregate_by_parcel = TRUE}),
#'   \code{decon_raw}, \code{decon_demixed}, \code{se_raw}, \code{se_demixed},
#'   \code{precision_weight}, \code{precision_weight_op},
#'   \code{precision_weight_raw}, \code{error_ac1}, \code{tau2_bio}, and the
#'   covariate columns for downstream \code{mixed_by()} modeling. Defaults to
#'   \code{TRUE} when \code{covariates = NULL} and \code{return_epochs = FALSE}.
#' @return A list with \code{coefficients} (per voxel/parcel, lag, and predictor
#'   term), \code{summary} (averaged across usable units per lag and term),
#'   optional \code{epochs}, optional \code{panel}, and \code{settings}.
#' @details
#'   \code{decon_epoch_regression()} bridges continuous regularized deconvolution
#'   and trial-level cognitive neuroscience models in five stages:
#'   \enumerate{
#'     \item \strong{Kernel Quadrature Readout:} At each requested lag \eqn{\tau},
#'       continuous deconvolved activity is integrated around each target event
#'       using \code{readout_kernel}. For \code{source_model = "duration"}, the
#'       readout window centers on the trial's midpoint
#'       \eqn{t_i + \tau + 0.5 d_i} and expands to \eqn{\max(w_{\mathrm{src}}, d_i)}.
#'     \item \strong{Non-Target Event De-Mixing (\code{demix_other_events}):}
#'       Using the frozen inverse operator \eqn{\mathbf{D}_r} and forward HRF
#'       design \eqn{\mathbf{G}_{r,0}}, the function subtracts resolution-predicted
#'       leakage from non-target events (e.g., cues or motor responses) using a
#'       condition-mean operator (\code{demix_ridge >= 0.03}) or a 2nd-order
#'       Tikhonov joint \eqn{(T, O)} operator (\code{demix_ridge < 0.03}).
#'     \item \strong{Target Cross-Talk & Duration De-Mixing (\code{demix_target_trials}):}
#'       Inverts the lag-specific target resolution matrix via a 2nd-order
#'       iterated Tikhonov operator (\code{solve_duration_demix_operator()}) with
#'       per-trial effective gain restoration to remove trial-to-trial overlap
#'       and divide out regularization-induced amplitude shrinkage. In
#'       \code{source_model = "duration"}, a centered-subspace penalty and
#'       rank-1 baseline adjustment keep constant baseline firing rates invariant
#'       to trial duration (preventing reaction-time confounding).
#'     \item \strong{Optional Multi-Lag Point-Spread De-Mixing (\code{demix_within_trial}):}
#'       When enabled across \eqn{\ge 3} lags, sharpens temporal specificity across
#'       lags using a 2nd-order Tikhonov lag-by-lag resolution inverse while
#'       preserving cross-trial operator covariance structure.
#'     \item \strong{Analytical Uncertainty, Precision Weights & Null-Mode-Damped ZCA Whitening:}
#'       Propagates observation-level AR(1) noise covariance through the composite
#'       linear readout/de-mixing operator \eqn{\mathbf{B}_{r,\tau}} to obtain exact
#'       per-trial standard errors (\code{se_demixed}), adjacent-trial error
#'       autocorrelation (\code{error_ac1}), Empirical-Bayes precision weights
#'       (\code{precision_weight}), covariance-propagated sandwich SEs
#'       (\code{se_weighted}), and null-mode-damped, Bartlett-adaptive within-run
#'       ZCA-whitened variables for \code{\link{mixed_by}}.
#'   }
#' @section Recommended Best-Practice Configuration:
#'   Across the 5-pillar MEDuSA Monte Carlo validation suite (amplitude recovery,
#'   temporal alignment, inter-region connectivity, factorial power surfaces, and
#'   biophysical Windkessel-Balloon + \code{fmrisim} scanner realism), the
#'   recommended end-to-end configuration is:
#'   \enumerate{
#'     \item \strong{Stage 0--1 (HRF Estimation):} Use
#'       \code{hrf_library(tr, preset = "medusa")} (0.5s peak-delay grid over
#'       4--8s with dispersions 0.9, 1.0, 1.1; 27 candidates) and
#'       \code{\link{estimate_crossrun_hrf}} with \code{noise_model = "ar1"},
#'       \code{evidence = "calibrated"}, \code{selection_mode = "soft_average"}
#'       (evidence-weighted HRF kernel averaging across locally supported
#'       shapes), and parcel labels in \code{spatial_groups}.
#'     \item \strong{Stage 2 (Deconvolution):} Use \code{\link{deconvolve_crossrun_hrf}}
#'       with \code{lambda = NULL} (adaptive prewhitened GCV tuning),
#'       \code{scaling = "noise_sd"} (AR(1) innovation noise SD rather than
#'       total \code{"zscore"} variance), and \code{penalty = "sobolev"}.
#'     \item \strong{Stage 3 (Epoch De-Mixing, Parcel Pre-Aggregation & Whitening):} Use
#'       \code{decon_epoch_regression(..., demix_other_events = TRUE,
#'       demix_target_trials = TRUE, restore_scale = TRUE, weighting = "gls",
#'       aggregate_by_parcel = TRUE, whiten_panel = TRUE, return_panel = TRUE)}.
#'       Enable \code{demix_within_trial = TRUE} with \code{lags} spaced at
#'       \eqn{1\times\mathrm{TR}} when estimating multi-lag waveforms \eqn{\beta(\tau)}
#'       (if \code{lags} is sampled finer than \eqn{1\times\mathrm{TR}},
#'       \code{within_trial_ridge} is automatically scaled by
#'       \eqn{(\mathrm{TR}/\Delta\tau)^2} to prevent sub-TR collinearity
#'       amplification), and set \code{source_model = "duration"}
#'       (\code{duration_scale = "rate"}) when modeling variable-duration (RT)
#'       epochs.
#'     \item \strong{Stage 4 (Group Mixed-Effects Inference):} Pass both the
#'       ZCA-whitened panel columns and precision weights to \code{\link{mixed_by}}:
#'       \code{mixed_by(panel, decon_demixed_whitened ~ x_whitened + (1 + x_whitened | id),
#'       weights = "precision_weight")}.
#'   }
#' @export
decon_epoch_regression <- function(
  deconvolution,
  covariates = NULL,
  events = NULL,
  target_condition = NULL,
  lags = seq(-2, 4, by = 1),
  readout_kernel = c("triangle", "linear", "boxcar", "gaussian"),
  readout_width = NULL,
  demix_other_events = TRUE,
  demix_target_trials = TRUE,
  demix_within_trial = FALSE,
  demix_ridge = 0.05,
  within_trial_ridge = 0.06,
  gain_floor = 0.05,
  restore_scale = FALSE,
  source_model = c("impulse", "window", "duration"),
  trial_durations = NULL,
  duration_scale = c("rate", "integral"),
  weighting = c("precision", "gls", "none"),
  tau2_mode = c("auto", "residual", "marginal"),
  aggregate_by_parcel = FALSE,
  parcel_groups = NULL,
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
  within_ridge_supplied <- !missing(within_trial_ridge)
  return_panel_supplied <- !missing(return_panel)
  readout_kernel <- match.arg(readout_kernel)
  source_model <- match.arg(source_model)
  duration_scale <- match.arg(duration_scale)
  weighting <- match.arg(weighting)
  tau2_mode <- match.arg(tau2_mode)
  if (!is.null(trial_durations) && !source_model_supplied) {
    source_model <- "duration"
  }
  if (!is.null(parcel_groups)) {
    aggregate_by_parcel <- TRUE
  }
  checkmate::assert_numeric(lags, min.len = 1L, any.missing = FALSE, finite = TRUE)
  checkmate::assert_string(target_condition, null.ok = TRUE)
  checkmate::assert_flag(demix_other_events)
  checkmate::assert_flag(demix_target_trials)
  checkmate::assert_flag(demix_within_trial)
  checkmate::assert_number(demix_ridge, lower = 1e-6, finite = TRUE)
  checkmate::assert_number(within_trial_ridge, lower = 1e-6, finite = TRUE)
  checkmate::assert_number(gain_floor, lower = 1e-4, finite = TRUE)
  checkmate::assert_flag(restore_scale)
  checkmate::assert_flag(aggregate_by_parcel)
  checkmate::assert_flag(whiten_panel)
  checkmate::assert_flag(return_epochs)
  checkmate::assert_flag(return_panel)
  if (is.null(covariates) && !return_panel_supplied && !isTRUE(return_epochs)) {
    return_panel <- TRUE
  }

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
  if (is.null(covariates)) {
    for (r in seq_len(n_runs)) {
      cov_list[[r]] <- matrix(numeric(0L), nrow = target_counts[r], ncol = 0L)
    }
    term_names <- character(0L)
    n_terms <- 0L
    eff_tau2_mode <- "marginal"
  } else if (is.list(covariates) && !is.data.frame(covariates)) {
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
    term_names <- colnames(cov_list[[1L]])
    n_terms <- length(term_names)
    eff_tau2_mode <- tau2_mode
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
    term_names <- colnames(cov_list[[1L]])
    n_terms <- length(term_names)
    eff_tau2_mode <- tau2_mode
  }

  n_voxels <- length(fit$voxel_ids)
  n_lags <- length(lags)
  tr <- fit$settings$tr

  # Resolve spatial parcel grouping and precompute cross-run inverse-variance voxel weights
  parcel_vox_list <- NULL
  parcel_weights_list <- NULL
  if (isTRUE(aggregate_by_parcel)) {
    if (!is.null(parcel_groups)) {
      if (!is.null(names(parcel_groups)) && all(fit$voxel_ids %in% names(parcel_groups))) {
        p_groups <- as.character(unname(parcel_groups[fit$voxel_ids]))
      } else {
        if (length(parcel_groups) != n_voxels) {
          stop("parcel_groups length (", length(parcel_groups),
            ") must equal number of voxels (", n_voxels, ")")
        }
        p_groups <- as.character(parcel_groups)
      }
    } else if (!is.null(fit$specification$spatial_groups)) {
      p_groups <- as.character(fit$specification$spatial_groups)
    } else {
      p_groups <- rep("parcel_1", n_voxels)
    }
    if (anyNA(p_groups) || any(!nzchar(p_groups))) {
      stop("parcel_groups must not contain NA or empty labels")
    }
    unit_ids <- unique(p_groups)
    n_units <- length(unit_ids)
    parcel_vox_list <- lapply(unit_ids, function( pid ) which(p_groups == pid))
    names(parcel_vox_list) <- unit_ids

    abs_s2_mat <- matrix(NA_real_, n_runs, n_voxels)
    rel_s2_mat <- matrix(NA_real_, n_runs, n_voxels)
    cv_r2_mat <- matrix(0.0, n_runs, n_voxels)
    for (r in seq_len(n_runs)) {
      sc_df_r <- deconvolution$scaling[[r]]
      ok_r <- sc_df_r$status == "ok"
      div_r0 <- sc_df_r$scale
      div_r0[!is.finite(div_r0) | div_r0 <= 0] <- 1
      nsd_r0 <- if (!is.null(sc_df_r$noise_sd)) sc_df_r$noise_sd else div_r0
      need_n0 <- ok_r & (!is.finite(nsd_r0) | nsd_r0 <= 0)
      if (any(need_n0)) nsd_r0[need_n0] <- div_r0[need_n0]
      rel_r0 <- nsd_r0 / div_r0
      rel_r0[!is.finite(rel_r0) | rel_r0 <= 0] <- 1
      abs_s2_mat[r, ok_r] <- nsd_r0[ok_r]^2
      rel_s2_mat[r, ok_r] <- rel_r0[ok_r]^2
      if (!is.null(fit$selections[[r]]$score)) {
        sc_v <- as.numeric(fit$selections[[r]]$score)
        sc_v[!is.finite(sc_v) | sc_v < 0] <- 0
        cv_r2_mat[r, ] <- sc_v
      }
    }
    parcel_weights_list <- vector("list", n_units)
    names(parcel_weights_list) <- unit_ids
    for (g in seq_len(n_units)) {
      vg <- parcel_vox_list[[g]]
      mean_abs_s2 <- colMeans(abs_s2_mat[, vg, drop = FALSE], na.rm = TRUE)
      mean_rel_s2 <- colMeans(rel_s2_mat[, vg, drop = FALSE], na.rm = TRUE)
      mean_cv_r2 <- colMeans(cv_r2_mat[, vg, drop = FALSE], na.rm = TRUE)
      mean_abs_s2[!is.finite(mean_abs_s2) | mean_abs_s2 <= 0] <- 1.0
      mean_rel_s2[!is.finite(mean_rel_s2) | mean_rel_s2 <= 0] <- 1.0
      med_abs <- stats::median(mean_abs_s2, na.rm = TRUE)
      med_rel <- stats::median(mean_rel_s2, na.rm = TRUE)
      s2_abs_fl <- pmax(mean_abs_s2, 0.10 * max(med_abs, 1e-8))
      s2_rel_fl <- pmax(mean_rel_s2, 0.10 * max(med_rel, 1e-8))
      med_r2 <- stats::median(mean_cv_r2, na.rm = TRUE)
      ev_mult <- if (any(mean_cv_r2 > 0)) {
        1.0 + mean_cv_r2 / (0.05 + max(med_r2, 1e-4))
      } else {
        rep(1.0, length(vg))
      }
      w_g <- if (isTRUE(restore_scale) || settings$scaling %in% c("demean", "none")) {
        ev_mult / s2_abs_fl
      } else {
        ev_mult / (s2_abs_fl * s2_rel_fl)
      }
      w_g[!is.finite(w_g) | w_g <= 0] <- 1.0
      parcel_weights_list[[g]] <- w_g / sum(w_g)
    }
  } else {
    unit_ids <- fit$voxel_ids
    n_units <- n_voxels
  }

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
  pw_epochs <- pw_op_epochs <- pw_raw_epochs <- error_ac1_epochs <- tau2_bio_epochs <- vector("list", n_lags)
  rel_epochs <- eb_y_epochs <- vector("list", n_lags)
  whitened_y_epochs <- se_whitened_epochs <- eb_whitened_y_epochs <- if (isTRUE(whiten_panel)) vector("list", n_lags) else NULL
  whitened_x_epochs <- if (isTRUE(whiten_panel) && n_terms > 0L) vector("list", n_lags) else NULL
  V_demix_blocks <- vector("list", n_lags)
  lag_names <- sprintf("lag_%g", lags)
  names(raw_epochs) <- names(demixed_epochs) <- names(se_raw_epochs) <- names(se_demixed_epochs) <- lag_names
  names(pw_epochs) <- names(pw_op_epochs) <- names(pw_raw_epochs) <- names(error_ac1_epochs) <- names(tau2_bio_epochs) <- lag_names
  names(rel_epochs) <- names(eb_y_epochs) <- lag_names
  if (isTRUE(whiten_panel)) {
    names(whitened_y_epochs) <- names(se_whitened_epochs) <- names(eb_whitened_y_epochs) <- lag_names
  }
  names(V_demix_blocks) <- lag_names
  for (li in seq_len(n_lags)) {
    na_mat <- matrix(NA_real_, total_trials, n_units, dimnames = list(NULL, unit_ids))
    raw_epochs[[li]] <- na_mat
    demixed_epochs[[li]] <- na_mat
    se_raw_epochs[[li]] <- na_mat
    se_demixed_epochs[[li]] <- na_mat
    pw_epochs[[li]] <- na_mat
    pw_op_epochs[[li]] <- na_mat
    pw_raw_epochs[[li]] <- na_mat
    error_ac1_epochs[[li]] <- na_mat
    tau2_bio_epochs[[li]] <- na_mat
    rel_epochs[[li]] <- na_mat
    eb_y_epochs[[li]] <- na_mat
    if (isTRUE(whiten_panel)) {
      whitened_y_epochs[[li]] <- na_mat
      se_whitened_epochs[[li]] <- na_mat
      eb_whitened_y_epochs[[li]] <- na_mat
      if (n_terms > 0L) {
        whitened_x_epochs[[li]] <- stats::setNames(
          replicate(n_terms, na_mat, simplify = FALSE),
          term_names
        )
      }
    }
    V_demix_blocks[[li]] <- stats::setNames(vector("list", n_runs), fit$run_ids)
    for (r in seq_len(n_runs)) V_demix_blocks[[li]][[r]] <- list()
  }
  gain_by_lag_voxel <- matrix(NA_real_, n_lags, n_units, dimnames = list(NULL, unit_ids))
  x_centered_list <- if (n_terms > 0L) {
    lapply(cov_list, function(x_r) sweep(x_r, 2L, colMeans(x_r, na.rm = TRUE), "-"))
  } else {
    cov_list
  }
  X_raw_pool <- do.call(rbind, cov_list)

  gain_accum <- matrix(0, n_lags, n_units)
  gain_count <- matrix(0L, n_lags, n_units)
  stratum_map <- vector("list", n_runs)
  sig_eff_mat <- matrix(NA_real_, n_runs, n_units, dimnames = list(fit$run_ids, unit_ids))
  center_obs <- !identical(settings$scaling, "none")
  lag_cov_accum <- stats::setNames(
    replicate(n_units, matrix(0, n_lags, n_lags), simplify = FALSE),
    unit_ids
  )
  lag_cov_trials <- integer(n_units)
  lag_res_accum <- stats::setNames(
    replicate(n_units, matrix(0, n_lags, n_lags), simplify = FALSE),
    unit_ids
  )
  lag_res_trials <- integer(n_units)

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

    strat <- hrf_stratum_kernels(fit, r)
    rho <- strat$rho
    keys <- strat$keys
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
    # Effective observation-noise SD in the working analysis units:
    # - If restore_scale = TRUE, deconvolved activity is multiplied back by div_r,
    #   so observation noise in restored units has SD = noise_sd_r.
    # - Otherwise, deconvolved activity is in standardized units (divided by div_r),
    #   so observation noise has SD = noise_sd_r / div_r (equals 1 under scaling = "noise_sd").
    sig_eff_r <- if (isTRUE(restore_scale)) noise_sd_r else (noise_sd_r / div_r)
    sig_eff_r[!is.finite(sig_eff_r) | sig_eff_r <= 0] <- 1

    if (isTRUE(restore_scale)) {
      act_r[, valid] <- sweep(act_r[, valid, drop = FALSE], 2L, div_r[valid], "*")
    }

    if (isTRUE(aggregate_by_parcel)) {
      act_parcel_r <- matrix(NA_real_, n, n_units, dimnames = list(NULL, unit_ids))
      valid_parcel <- logical(n_units)
      sig_eff_parcel_r <- rep(1.0, n_units)
      rho_parcel <- rep(0.0, n_units)
      keys_parcel <- character(n_units)
      kernels_parcel <- vector("list", n_units)
      ar_parcel <- vector("list", n_units)

      for (g in seq_len(n_units)) {
        vox_g <- parcel_vox_list[[g]]
        ok_g <- valid[vox_g]
        if (!any(ok_g)) next
        v_ok <- vox_g[ok_g]
        valid_parcel[g] <- TRUE
        u_g <- parcel_weights_list[[g]][ok_g]
        u_g <- u_g / sum(u_g)
        act_parcel_r[, g] <- as.numeric(act_r[, v_ok, drop = FALSE] %*% u_g)

        seff_v <- sig_eff_r[v_ok]
        s2_indep <- sum((u_g^2) * (seff_v^2))
        s2_ceil <- (sum(u_g * seff_v))^2
        if (length(v_ok) == 1L || length(rows) < 4L) {
          sig_eff_parcel_r[g] <- sqrt(s2_indep)
        } else {
          D2_vox <- diff(act_r[rows, v_ok, drop = FALSE], differences = 2L)
          d2_p <- as.numeric(D2_vox %*% u_g)
          mad_d2_v <- apply(D2_vox, 2L, stats::mad, na.rm = TRUE)
          mad_d2_p <- stats::mad(d2_p, na.rm = TRUE)
          denom_d2 <- sum(u_g * mad_d2_v)
          if (is.finite(mad_d2_p) && is.finite(denom_d2) && denom_d2 > 1e-12) {
            r_emp <- (mad_d2_p / denom_d2)^2
            s2_eff_g <- min(s2_ceil, max(s2_indep, r_emp * s2_ceil))
          } else {
            s2_eff_g <- s2_indep
          }
          sig_eff_parcel_r[g] <- sqrt(max(s2_eff_g, 1e-12))
        }

        u_keys <- unique(keys[v_ok])
        if (length(u_keys) == 1L) {
          keys_parcel[g] <- u_keys[1L]
          rho_parcel[g] <- rho[v_ok[1L]]
          kernels_parcel[[g]] <- strat$get_kernel(v_ok[1L])
          ar_parcel[[g]] <- if (!is.null(strat$get_ar)) strat$get_ar(v_ok[1L]) else rho[v_ok[1L]]
        } else {
          rho_g <- sum(u_g * rho[v_ok])
          psi1_g <- if (!is.null(strat$psi1)) sum(u_g * strat$psi1[v_ok]) else 0.0
          psi2_g <- if (!is.null(strat$psi2)) sum(u_g * strat$psi2[v_ok]) else 0.0
          k_mat <- vapply(v_ok, function(vi) as.numeric(strat$get_kernel(vi)), numeric(nrow(fit$kernels)))
          k_g <- as.numeric(k_mat %*% u_g)
          k_sum <- sum(k_g)
          if (is.finite(k_sum) && abs(k_sum) > 1e-6) k_g <- k_g / k_sum
          rho_parcel[g] <- rho_g
          kernels_parcel[[g]] <- matrix(k_g, ncol = 1L)
          keys_parcel[g] <- paste0("parcel_", unit_ids[g], ":", round(rho_g, 6))
          ar_parcel[[g]] <- if (abs(psi1_g) > 1e-10 || abs(psi2_g) > 1e-10) c(rho_g, psi1_g, psi2_g) else rho_g
        }
      }
      act_r <- act_parcel_r
      valid <- valid_parcel
      sig_eff_r <- sig_eff_parcel_r
      rho <- rho_parcel
      keys <- keys_parcel
      get_unit_kernel <- function(idx) kernels_parcel[[idx]]
      get_unit_ar <- function(idx) ar_parcel[[idx]]
    } else {
      get_unit_kernel <- strat$get_kernel
      get_unit_ar <- if (!is.null(strat$get_ar)) strat$get_ar else function(idx) rho[idx]
    }
    stratum_map[[r]] <- keys
    sig_eff_mat[r, ] <- sig_eff_r

    # Build readout matrices:
    # - S_lags[[li]]: standard readout centered at onset + lag (used for raw readouts
    #   and for impulse/window models).
    # - S_lags_model[[li]]: for source_model = "duration", centers the readout at the
    #   trial boxcar midpoint (onset + lag + 0.5 * d_i) with trial-specific width
    #   max(w_src, d_i), aligning sustained processes of varying duration at lag = 0.
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
      ar <- get_unit_ar(voxels[1L])
      hrf_k <- get_unit_kernel(voxels[1L])
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

      L_joint <- if (demix_other_events && length(o_idx) > 0L) {
        build_other_demix_operator(S_all_0, DG_all_0, ev, o_idx, t_idx, demix_ridge,
          source_model = source_model, demix_target_trials = demix_target_trials)
      } else {
        NULL
      }
      L_other <- if (!is.null(L_joint)) L_joint$L_other else NULL
      L_target_0 <- if (!is.null(L_joint)) L_joint$L_target_0 else NULL
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
          if (!is.null(L_other)) {
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

      ok_tau_list <- vector("list", n_lags)
      B_demix_list <- vector("list", n_lags)
      DG_tau_list <- vector("list", n_lags)

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
        ok_tau_list[[li]] <- ok_tau
        if (!any(ok_tau)) next
        m_ok <- sum(ok_tau)
        Y_clean <- demixed_epochs[[li]][row_span, voxels, drop = FALSE]
        S_clean_tau <- S_tau[ok_tau, , drop = FALSE]

        max_run_time <- offset + (n - 1) * tr
        if (source_model == "duration") {
          shifted_onsets <- pmin(max_run_time, pmax(offset, ev$onset[t_idx] + lags[li]))
          G_tau_unc <- hrf_event_design_cpp(n, tr, offset, shifted_onsets, dur_target,
            seq_len(m_target), m_target, hrf_k)[[1L]][rows, , drop = FALSE]
          if (duration_scale == "integral" && any(dur_target > 0)) {
            d_pos <- dur_target > 0
            G_tau_unc[, d_pos] <- sweep(G_tau_unc[, d_pos, drop = FALSE], 2L, dur_target[d_pos] / tr, "/")
          }
        } else {
          # Local target resolution matrix K_tau for unit sources at target onset + lag
          dur_tau <- if (source_model == "window") rep(w_src, m_target) else rep(0, m_target)
          shift_center <- if (source_model == "window") lags[li] - 0.5 * w_src else lags[li]
          shifted_onsets <- pmin(max_run_time, pmax(offset, ev$onset[t_idx] + shift_center))
          G_tau_unc <- hrf_event_design_cpp(n, tr, offset, shifted_onsets, dur_tau,
            seq_len(m_target), m_target, hrf_k)[[1L]][rows, , drop = FALSE]
          if (source_model == "window" && w_src > 0) {
            G_tau_unc <- G_tau_unc / (w_src / tr)
          }
        }
        G_tau <- if (settings$scaling != "none") {
          sweep(G_tau_unc, 2L, colMeans(G_tau_unc), "-")
        } else {
          G_tau_unc
        }
        DG_tau_list[[li]] <- operator %*% G_tau

        if (!is.null(L_other)) {
          R_other_tau <- S_tau[ok_tau, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
          Y_clean[ok_tau, ] <- Y_clean[ok_tau, , drop = FALSE] - R_other_tau %*% a_other_hat
          S_clean_tau <- S_clean_tau - R_other_tau %*% L_other
        }
        B_clean_tau <- S_clean_tau %*% operator

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
          if (identical(source_model, "impulse") && !isTRUE(demix_within_trial) &&
              abs(lags[li]) <= 1e-6 && !is.null(L_target_0)) {
            L_t0_ok <- L_target_0[ok_tau, , drop = FALSE]
            Y_clean[ok_tau, ] <- L_t0_ok %*% act_r[, voxels, drop = FALSE]
            B_demix_tau <- L_t0_ok %*% operator
          } else {
            Y_sub <- Y_clean[ok_tau, , drop = FALSE]
            B_sub_tau <- B_clean_tau
            if (source_model == "duration" && !is.null(alpha0_hat) && !is.null(DG_target_0_sum) && abs(lags[li]) > 1e-8) {
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
          }
        } else {
          B_demix_tau <- B_clean_tau
        }
        demixed_epochs[[li]][row_span, voxels] <- Y_clean
        B_demix_list[[li]] <- B_demix_tau
      }

      # Cross-lag point-spread resolution matrix C_sym and optional within-trial de-mixing
      ok_common <- NULL
      m_common <- 0L
      if (n_lags >= 2L) {
        has_all_lags <- all(vapply(ok_tau_list, function(x) !is.null(x) && length(x) == m_target, logical(1L)))
        if (has_all_lags) {
          ok_common <- Reduce("&", ok_tau_list)
          m_common <- sum(ok_common)
          if (m_common >= 2L) {
            K_lag <- matrix(0, n_lags, n_lags)
            for (l1 in seq_len(n_lags)) {
              S_l1_sub <- S_lags_model[[l1]][ok_common, , drop = FALSE]
              for (l2 in seq_len(n_lags)) {
                DG_l2_sub <- DG_tau_list[[l2]][, ok_common, drop = FALSE]
                k_diag <- rowSums(S_l1_sub * t(DG_l2_sub))
                K_lag[l1, l2] <- mean(k_diag[is.finite(k_diag)])
              }
            }
            d_lag <- pmax(diag(K_lag), gain_floor)
            D_inv_sqrt <- diag(1 / sqrt(d_lag), n_lags)
            C_lag <- D_inv_sqrt %*% K_lag %*% D_inv_sqrt
            C_sym <- 0.5 * (C_lag + t(C_lag))
            for (vi in voxels) {
              lag_res_accum[[vi]] <- lag_res_accum[[vi]] + m_common * C_sym
              lag_res_trials[vi] <- lag_res_trials[vi] + m_common
            }

            if (isTRUE(demix_within_trial) && n_lags >= 3L) {
              dt_lag <- if (n_lags > 1L) min(diff(sort(unique(lags)))) else tr
              dt_ref <- max(tr, 1.0)
              sub_tr_ratio <- if (is.finite(dt_lag) && dt_lag > 0 && dt_lag < dt_ref) {
                (dt_ref / max(dt_lag, 0.1 * dt_ref))^2
              } else {
                1.0
              }
              eff_wt_ridge <- if (!within_ridge_supplied && sub_tr_ratio > 1.0) {
                within_trial_ridge * sub_tr_ratio
              } else {
                within_trial_ridge
              }
              neumann_wt <- 0.6 / sub_tr_ratio
              W_1 <- solve(crossprod(C_sym) + eff_wt_ridge * diag(n_lags), t(C_sym))
              W_sym <- W_1 + neumann_wt * W_1 %*% (diag(n_lags) - C_sym %*% W_1)
              R_mat <- W_sym %*% C_sym
              diag_base <- pmax(diag(R_mat), gain_floor)
              W_lag_unit <- sweep(W_sym, 1L, diag_base, "/")
              if (is.finite(dt_lag) && dt_lag > 0 && dt_lag < dt_ref) {
                # Under sub-second lag spacing, single-lag recovery normalizes each lag by
                # its 1-bin impulse gain while smooth neural waveforms span a continuous
                # reference lobe (~0.72s). Scale W_lag by the interior reference-lobe gain
                # ratio so multi-lag de-mixing sharpens resolution without inflating waveform
                # amplitude or distorting boundary-lag correlation.
                sigma_lobe <- 0.72 * sqrt(max(0.0, 1.0 - (dt_lag / dt_ref)^2))
                lag_dist <- abs(outer(lags, lags, "-"))
                L_lobe <- exp(-0.5 * (lag_dist / max(sigma_lobe, 1e-6))^2)
                lobe_ratios <- rowSums(R_mat * L_lobe) / diag_base
                mid_idx <- if (n_lags >= 3L) seq.int(2L, n_lags - 1L) else seq_len(n_lags)
                lobe_scale <- max(1.0, stats::median(lobe_ratios[mid_idx], na.rm = TRUE))
                diag_gain <- diag_base * lobe_scale
              } else {
                diag_gain <- diag_base
              }
              W_lag <- sweep(W_sym, 1L, diag_gain, "/")

              idx_common <- which(ok_common)
              Y_common_old <- lapply(seq_len(n_lags), function(li) {
                demixed_epochs[[li]][row_span[idx_common], voxels, drop = FALSE]
              })
              B_common_old <- lapply(seq_len(n_lags), function(li) {
                ok_sub_li <- ok_common[ok_tau_list[[li]]]
                B_demix_list[[li]][ok_sub_li, , drop = FALSE]
              })
              n_obs_b <- ncol(B_common_old[[1L]])
              for (l1 in seq_len(n_lags)) {
                Y_new_l1 <- matrix(0, m_common, length(voxels))
                B_new_l1 <- matrix(0, m_common, n_obs_b)
                for (l2 in seq_len(n_lags)) {
                  w_12 <- W_lag[l1, l2]
                  if (abs(w_12) > 0) {
                    Y_new_l1 <- Y_new_l1 + w_12 * Y_common_old[[l2]]
                    B_new_l1 <- B_new_l1 + w_12 * B_common_old[[l2]]
                  }
                }
                demixed_epochs[[l1]][row_span[idx_common], voxels] <- Y_new_l1
                ok_sub_l1 <- ok_common[ok_tau_list[[l1]]]
                B_demix_list[[l1]][ok_sub_l1, ] <- B_new_l1
              }
            }
          }
        }
      }

      U_demix_list <- vector("list", n_lags)
      for (li in seq_len(n_lags)) {
        ok_tau <- ok_tau_list[[li]]
        if (is.null(ok_tau) || !any(ok_tau)) next
        m_ok <- sum(ok_tau)
        B_demix_tau <- B_demix_list[[li]]

        # Compute dimensionless trial innovation factor and error covariance for stratum (r, key, li)
        U_demix <- compute_ar1_operator_factor(B_demix_tau, rows, ar, center_obs = center_obs)
        U_demix_list[[li]] <- U_demix
        V_demix <- crossprod(U_demix)
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

      # Exact L x L cross-lag trial error covariance across lags in a single BLAS crossprod
      if (n_lags >= 2L && m_common >= 2L && !any(vapply(U_demix_list, is.null, logical(1L)))) {
        n_obs_u <- nrow(U_demix_list[[1L]])
        U_lag_mat <- matrix(0, n_obs_u * m_common, n_lags)
        for (li in seq_len(n_lags)) {
          ok_sub_li <- ok_common[ok_tau_list[[li]]]
          U_sub <- U_demix_list[[li]][, ok_sub_li, drop = FALSE]
          U_cent <- sweep(U_sub, 1L, rowMeans(U_sub), "-")
          U_lag_mat[, li] <- as.numeric(U_cent)
        }
        V_lag_unit <- crossprod(U_lag_mat) / max(1L, m_common - 1L)
        for (vi in voxels) {
          lag_cov_accum[[vi]] <- lag_cov_accum[[vi]] + m_common * (sig_eff_r[vi]^2) * V_lag_unit
          lag_cov_trials[vi] <- lag_cov_trials[vi] + m_common
        }
      }
    }
  }

  has_gain <- gain_count > 0L
  gain_by_lag_voxel[has_gain] <- gain_accum[has_gain] / gain_count[has_gain]

  # Center covariates within run for pooled regression
  X_pool <- do.call(rbind, x_centered_list)
  coef_rows <- vector("list", n_lags * n_units * max(1L, n_terms))
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

    for (vi in seq_len(n_units)) {
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

      ok_v <- if (n_terms > 0L) {
        is.finite(y_raw_v) & is.finite(y_dem_v) & (rowSums(!is.finite(X_pool)) == 0L)
      } else {
        is.finite(y_raw_v) & is.finite(y_dem_v)
      }
      n_ok <- sum(ok_v)
      mean_se_v <- if (any(ok_s_dem)) mean(s_dem_v[ok_s_dem]) else NA_real_
      ac1_v <- error_ac1_epochs[[li]][, vi]
      mean_ac1_v <- if (any(is.finite(ac1_v))) mean(ac1_v[is.finite(ac1_v) & ac1_v != 0]) else NA_real_
      if (!is.finite(mean_ac1_v)) mean_ac1_v <- 0

      if (n_ok <= n_terms + 2L) {
        if (any(ok_s_dem)) {
          pw_epochs[[li]][ok_s_dem, vi] <- pw_op_epochs[[li]][ok_s_dem, vi]
        }
        if (n_terms > 0L) {
          for (ti in seq_len(n_terms)) {
            row_k <- row_k + 1L
            df_row <- data.frame(
              voxel = unit_ids[vi], lag = lags[li], term = term_names[ti],
              beta_raw = NA_real_, beta_corrected = NA_real_,
              beta_unweighted = NA_real_, beta_weighted = NA_real_, beta_gls = NA_real_,
              se_corrected = NA_real_, se_unweighted = NA_real_, se_weighted = NA_real_, se_gls = NA_real_,
              r_raw = NA_real_, r_demixed = NA_real_,
              t_stat = NA_real_, p_value = NA_real_,
              t_unweighted = NA_real_, t_weighted = NA_real_, t_gls = NA_real_,
              df = n_ok - n_terms,
              tau2_bio = NA_real_,
              mean_gain = gain_by_lag_voxel[li, vi],
              mean_se_trial = mean_se_v,
              mean_error_ac1 = mean_ac1_v,
              stringsAsFactors = FALSE
            )
            if (isTRUE(aggregate_by_parcel)) df_row$parcel <- unit_ids[vi]
            coef_rows[[row_k]] <- df_row
          }
        } else {
          row_k <- row_k + 1L
          df_row <- data.frame(
            voxel = unit_ids[vi], lag = lags[li], term = NA_character_,
            beta_raw = NA_real_, beta_corrected = NA_real_,
            beta_unweighted = NA_real_, beta_weighted = NA_real_, beta_gls = NA_real_,
            se_corrected = NA_real_, se_unweighted = NA_real_, se_weighted = NA_real_, se_gls = NA_real_,
            r_raw = NA_real_, r_demixed = NA_real_,
            t_stat = NA_real_, p_value = NA_real_,
            t_unweighted = NA_real_, t_weighted = NA_real_, t_gls = NA_real_,
            df = max(0L, n_ok - n_runs),
            tau2_bio = NA_real_,
            mean_gain = gain_by_lag_voxel[li, vi],
            mean_se_trial = mean_se_v,
            mean_error_ac1 = mean_ac1_v,
            stringsAsFactors = FALSE
          )
          if (isTRUE(aggregate_by_parcel)) df_row$parcel <- unit_ids[vi]
          coef_rows[[row_k]] <- df_row
        }
        next
      }

      # Assemble per-run operator covariance blocks Omega_r = sigma_{e,r,v}^2 * V_{demix,cent}
      # and compute both covariate-agnostic marginal and OLS-residual biological variance
      # estimators tau2_bio. Combining tau2_bio with the diagonal of Omega_r yields
      # Empirical-Bayes trial precision weights w_i = 1 / (tau2_eff + s_{demix,i}^2),
      # preventing over-weighting high-gain trials when biological variability dominates.
      Omega_r_list <- vector("list", n_runs)
      rs_ok_list <- vector("list", n_runs)
      xOx_sum <- if (n_terms > 0L) matrix(0, n_terms, n_terms) else NULL
      tr_Omega <- 0
      ss_marg <- 0
      tr_Omega_marg <- 0
      df_marg <- 0L
      ss_diff <- 0
      exp_diff_Om <- 0
      n_diff <- 0L

      for (r in seq_len(n_runs)) {
        rs <- starts[r]:ends[r]
        ok_r_local <- which(ok_v[rs])
        m_r <- length(ok_r_local)
        if (!m_r) next
        rs_r <- rs[ok_r_local]
        rs_ok_list[[r]] <- rs_r
        key_rv <- if (!is.null(stratum_map[[r]])) stratum_map[[r]][vi] else NULL
        blk_r <- if (!is.null(key_rv)) V_demix_blocks[[li]][[r]][[key_rv]] else NULL
        if (!is.null(blk_r)) {
          ok_tau_idx <- which(blk_r$ok_tau)
          pos_in_blk <- match(ok_r_local, ok_tau_idx)
          valid_pos <- !is.na(pos_in_blk)
          if (all(valid_pos)) {
            Omega_r <- (sig_eff_mat[r, vi]^2) * blk_r$V_cent[pos_in_blk, pos_in_blk, drop = FALSE]
          } else {
            s2_sub <- s2_dem_v[rs_r]
            s2_sub[!is.finite(s2_sub)] <- 1
            Omega_r <- diag(s2_sub, m_r)
          }
        } else {
          s2_sub <- s2_dem_v[rs_r]
          s2_sub[!is.finite(s2_sub)] <- 1
          Omega_r <- diag(s2_sub, m_r)
        }
        Omega_r_list[[r]] <- Omega_r
        tr_r <- sum(diag(Omega_r))
        tr_Omega <- tr_Omega + tr_r
        if (m_r >= 2L) {
          yr_dem_r <- demixed_epochs[[li]][rs_r, vi]
          yr_c_r <- yr_dem_r - mean(yr_dem_r)
          ss_marg <- ss_marg + sum(yr_c_r^2)
          tr_Omega_marg <- tr_Omega_marg + ((m_r - 1.0) / m_r) * tr_r
          df_marg <- df_marg + (m_r - 1L)
          dy_r <- diff(yr_dem_r)
          ss_diff <- ss_diff + sum(dy_r^2)
          diag_Om <- diag(Omega_r)
          off1_Om <- Omega_r[cbind(2:m_r, 1:(m_r - 1L))]
          exp_diff_Om <- exp_diff_Om + sum(diag_Om[-1L] + diag_Om[-m_r] - 2.0 * off1_Om)
          n_diff <- n_diff + (m_r - 1L)
        }
        if (n_terms > 0L) {
          Xr_ok <- X_pool[rs_r, , drop = FALSE]
          xOx_sum <- xOx_sum + crossprod(Xr_ok, Omega_r %*% Xr_ok)
        }
      }

      s2_ok <- s2_dem_v[ok_v]
      s2_med <- if (any(is.finite(s2_ok) & s2_ok > 0)) stats::median(s2_ok[is.finite(s2_ok) & s2_ok > 0]) else 1
      tau2_var <- if (df_marg >= 1L) max(0, (ss_marg - tr_Omega_marg) / df_marg) else 0
      tau2_diff <- if (n_diff >= 1L) max(0, (ss_diff - exp_diff_Om) / (2.0 * n_diff)) else 0
      mad_y <- stats::mad(y_dem_v[ok_v], constant = 1.4826, na.rm = TRUE)
      tau2_mad <- if (is.finite(mad_y)) max(0, mad_y^2 - s2_med) else 0
      tau2_marginal <- 0.45 * tau2_var + 0.35 * tau2_diff + 0.20 * tau2_mad

      if (n_terms > 0L) {
        X_ok <- X_pool[ok_v, , drop = FALSE]
        xtx_raw <- crossprod(X_ok) + 1e-8 * diag(n_terms)
        xtx_inv_raw <- solve(xtx_raw)
        beta_raw_v <- as.numeric(xtx_inv_raw %*% crossprod(X_ok, y_raw_v[ok_v]))
        beta_unw_v <- as.numeric(xtx_inv_raw %*% crossprod(X_ok, y_dem_v[ok_v]))

        resid_unw <- y_dem_v[ok_v] - as.numeric(X_ok %*% beta_unw_v)
        df_v <- max(1L, n_ok - n_terms - n_runs)
        sigma2_unw <- sum(resid_unw^2) / df_v
        se_unw_naive <- sqrt(pmax(diag(xtx_inv_raw) * sigma2_unw, .Machine$double.eps))

        tr_proj_Omega <- max(0, tr_Omega - sum(diag(xtx_inv_raw %*% xOx_sum)))
        tau2_resid <- max(0, (sum(resid_unw^2) - tr_proj_Omega) / df_v)
        tau2_bio <- if (identical(eff_tau2_mode, "marginal")) {
          tau2_marginal
        } else if (identical(eff_tau2_mode, "residual")) {
          tau2_resid
        } else {
          0.65 * tau2_resid + 0.35 * tau2_marginal
        }
      } else {
        df_v <- max(1L, n_ok - n_runs)
        tau2_bio <- tau2_marginal
      }
      tau2_bio_epochs[[li]][, vi] <- tau2_bio
      tau2_eff <- max(tau2_bio, 0.10 * s2_med)

      if (any(ok_s_dem)) {
        w_eb_all <- 1 / (tau2_eff + s2_dem_v[ok_s_dem])
        pw_epochs[[li]][ok_s_dem, vi] <- w_eb_all / mean(w_eb_all)
        rel_epochs[[li]][ok_s_dem, vi] <- tau2_eff / (tau2_eff + s2_dem_v[ok_s_dem])
      } else {
        pw_epochs[[li]][ok_v, vi] <- 1
        rel_epochs[[li]][ok_v, vi] <- 1
      }
      w_v_all <- pw_epochs[[li]][, vi]
      w_v_all[!is.finite(w_v_all) | w_v_all <= 0] <- 1
      rel_v_all <- rel_epochs[[li]][, vi]
      rel_v_all[!is.finite(rel_v_all) | rel_v_all < 0] <- 1
      for (r in seq_len(n_runs)) {
        rs_ok <- rs_ok_list[[r]]
        if (!length(rs_ok)) next
        yr_raw <- demixed_epochs[[li]][rs_ok, vi]
        mu_r <- mean(yr_raw)
        eb_y_epochs[[li]][rs_ok, vi] <- mu_r + rel_v_all[rs_ok] * (yr_raw - mu_r)
      }

      if (n_terms == 0L) {
        if (isTRUE(whiten_panel)) {
          for (r in seq_len(n_runs)) {
            rs_ok <- rs_ok_list[[r]]
            m_r <- length(rs_ok)
            if (!m_r) next
            yr_raw <- demixed_epochs[[li]][rs_ok, vi]
            yr_cent <- yr_raw - mean(yr_raw)
            Omega_r <- Omega_r_list[[r]]
            tau2_zca <- max(tau2_eff, 0.15 * s2_med)
            Sigma_tot_r <- Omega_r + diag(tau2_zca, m_r)
            Phi_r <- compute_zca_whitening_matrix(Sigma_tot_r, eig_floor = 0.05, offdiag_scale = 0.85)
            yr_zca <- as.numeric(Phi_r %*% yr_cent)
            s2_zca_r <- pmax(0, rowSums((Phi_r %*% Omega_r) * Phi_r))
            se_zca_r <- sqrt(s2_zca_r)
            rel_zca_r <- tau2_eff / (tau2_eff + s2_zca_r)
            whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + yr_zca
            se_whitened_epochs[[li]][rs_ok, vi] <- se_zca_r
            eb_whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + rel_zca_r * yr_zca
          }
        }
        row_k <- row_k + 1L
        df_row <- data.frame(
          voxel = unit_ids[vi], lag = lags[li], term = NA_character_,
          beta_raw = NA_real_, beta_corrected = NA_real_,
          beta_unweighted = NA_real_, beta_weighted = NA_real_, beta_gls = NA_real_,
          se_corrected = NA_real_, se_unweighted = NA_real_, se_weighted = NA_real_, se_gls = NA_real_,
          r_raw = NA_real_, r_demixed = NA_real_,
          t_stat = NA_real_, p_value = NA_real_,
          t_unweighted = NA_real_, t_weighted = NA_real_, t_gls = NA_real_,
          df = df_v,
          tau2_bio = tau2_bio,
          mean_gain = gain_by_lag_voxel[li, vi],
          mean_se_trial = mean_se_v,
          mean_error_ac1 = mean_ac1_v,
          stringsAsFactors = FALSE
        )
        if (isTRUE(aggregate_by_parcel)) df_row$parcel <- unit_ids[vi]
        coef_rows[[row_k]] <- df_row
        next
      }

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
        tau2_zca <- max(tau2_eff, 0.15 * s2_med)
        Sigma_tot_r <- Omega_r + diag(tau2_zca, m_r)
        Sigma_tot_list[[r]] <- list(idx = idx_span, Sigma = Sigma_tot_r, wr = wr)

        yr_cent <- yr_raw - mean(yr_raw)
        Xr_cent <- sweep(Xr_raw, 2L, colMeans(Xr_raw), "-")
        rho_x_r <- if (m_r >= 4L) {
          max(vapply(seq_len(n_terms), function(ti) {
            xv <- Xr_cent[, ti]
            if (stats::sd(xv) > 1e-10) abs(stats::cor(xv[-1L], xv[-m_r])) else 0.0
          }, numeric(1L)), na.rm = TRUE)
        } else {
          1.0
        }
        if (!is.finite(rho_x_r)) rho_x_r <- 0.5
        gamma_r <- min(1.0, max(0.25, 1.8 * rho_x_r))
        Phi_r <- compute_zca_whitening_matrix(Sigma_tot_r, eig_floor = 0.05, offdiag_scale = gamma_r)
        yr_zca <- as.numeric(Phi_r %*% yr_cent)
        Xr_zca <- Phi_r %*% Xr_cent
        if (isTRUE(whiten_panel)) {
          s2_zca_r <- pmax(0, rowSums((Phi_r %*% Omega_r) * Phi_r))
          se_zca_r <- sqrt(s2_zca_r)
          rel_zca_r <- tau2_eff / (tau2_eff + s2_zca_r)
          whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + yr_zca
          se_whitened_epochs[[li]][rs_ok, vi] <- se_zca_r
          eb_whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + rel_zca_r * yr_zca
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
        df_row <- data.frame(
          voxel = unit_ids[vi], lag = lags[li], term = term_names[ti],
          beta_raw = beta_raw_v[ti], beta_corrected = beta_corr_v[ti],
          beta_unweighted = beta_unw_v[ti], beta_weighted = beta_pw_v[ti], beta_gls = beta_gls_v[ti],
          se_corrected = se_corr_v[ti], se_unweighted = se_unw_naive[ti],
          se_weighted = se_pw_v[ti], se_gls = se_gls_v[ti],
          r_raw = r_raw_val, r_demixed = r_dem_val,
          t_stat = t_corr_vec[ti], p_value = p_corr_vec[ti],
          t_unweighted = t_unw_vec[ti], t_weighted = t_pw_vec[ti], t_gls = t_gls_vec[ti],
          df = df_v,
          tau2_bio = tau2_bio,
          mean_gain = gain_by_lag_voxel[li, vi],
          mean_se_trial = mean_se_v,
          mean_error_ac1 = mean_ac1_v,
          stringsAsFactors = FALSE
        )
        if (isTRUE(aggregate_by_parcel)) df_row$parcel <- unit_ids[vi]
        coef_rows[[row_k]] <- df_row
      }
    }
  }

  coef_df <- do.call(rbind, coef_rows[seq_len(row_k)])

  lag_names <- as.character(lags)
  lag_error_cov <- stats::setNames(vector("list", n_units), unit_ids)
  lag_error_cor <- stats::setNames(vector("list", n_units), unit_ids)
  lag_resolution <- stats::setNames(vector("list", n_units), unit_ids)
  for (vi in seq_len(n_units)) {
    if (n_lags >= 2L && lag_cov_trials[vi] > 0L) {
      V_u <- lag_cov_accum[[vi]] / lag_cov_trials[vi]
      V_u <- 0.5 * (V_u + t(V_u))
      d_u <- sqrt(pmax(diag(V_u), .Machine$double.eps))
      R_u <- V_u / outer(d_u, d_u)
      R_u[!is.finite(R_u)] <- 0
      diag(R_u) <- 1
    } else {
      V_u <- diag(1, n_lags)
      R_u <- diag(1, n_lags)
    }
    if (n_lags >= 2L && lag_res_trials[vi] > 0L) {
      C_u <- lag_res_accum[[vi]] / lag_res_trials[vi]
      C_u <- 0.5 * (C_u + t(C_u))
    } else {
      C_u <- diag(1, n_lags)
    }
    dimnames(V_u) <- dimnames(R_u) <- dimnames(C_u) <- list(lag_names, lag_names)
    lag_error_cov[[vi]] <- V_u
    lag_error_cor[[vi]] <- R_u
    lag_resolution[[vi]] <- C_u
  }

  if (nrow(coef_df) > 0L) {
    coef_df$beta_multilag <- coef_df$beta_corrected
    coef_df$se_multilag <- coef_df$se_corrected
    if (n_lags >= 3L && n_terms > 0L) {
      for (vi in seq_len(n_units)) {
        u_id <- unit_ids[vi]
        C_u <- lag_resolution[[vi]]
        R_u <- lag_error_cor[[vi]]
        for (ti in seq_len(n_terms)) {
          idx_vt <- which(coef_df$voxel == u_id & coef_df$term == term_names[ti])
          if (length(idx_vt) == n_lags) {
            ord_l <- match(lags, coef_df$lag[idx_vt])
            idx_vt <- idx_vt[ord_l]
            b_v <- coef_df$beta_corrected[idx_vt]
            se_v <- coef_df$se_corrected[idx_vt]
            if (all(is.finite(b_v)) && all(is.finite(se_v) & se_v > 0)) {
              Sigma_b <- outer(se_v, se_v) * R_u
              sm_res <- medusa_smooth_deblur_trajectory(
                beta_vec = b_v,
                Sigma_beta = Sigma_b,
                lags = lags,
                C_sym = C_u,
                tr = fit$settings$tr,
                demix_within_trial = isTRUE(demix_within_trial) || !isTRUE(demix_target_trials),
                within_trial_ridge = within_trial_ridge,
                trajectory_mode = "auto"
              )
              coef_df$beta_multilag[idx_vt] <- sm_res$beta
              coef_df$se_multilag[idx_vt] <- sm_res$se
            }
          }
        }
      }
    }
    attr(coef_df, "lag_error_cov") <- lag_error_cov
    attr(coef_df, "lag_error_cor") <- lag_error_cor
    attr(coef_df, "lag_resolution") <- lag_resolution
  }

  fmean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  summary_rows <- list()
  for (li in seq_len(n_lags)) {
    if (n_terms > 0L) {
      for (ti in seq_len(n_terms)) {
        sub_df <- coef_df[coef_df$lag == lags[li] & coef_df$term == term_names[ti], , drop = FALSE]
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          lag = lags[li], term = term_names[ti],
          mean_beta_raw = fmean(sub_df$beta_raw),
          mean_beta_corrected = fmean(sub_df$beta_corrected),
          mean_beta_multilag = fmean(sub_df$beta_multilag),
          mean_beta_unweighted = fmean(sub_df$beta_unweighted),
          mean_beta_weighted = fmean(sub_df$beta_weighted),
          mean_beta_gls = fmean(sub_df$beta_gls),
          mean_r_raw = fmean(sub_df$r_raw),
          mean_r_demixed = fmean(sub_df$r_demixed),
          mean_tau2_bio = fmean(sub_df$tau2_bio),
          mean_gain = fmean(sub_df$mean_gain),
          mean_se_trial = fmean(sub_df$mean_se_trial),
          mean_error_ac1 = fmean(sub_df$mean_error_ac1),
          n_voxels = sum(is.finite(sub_df$beta_corrected)),
          stringsAsFactors = FALSE
        )
      }
    } else {
      sub_df <- coef_df[coef_df$lag == lags[li], , drop = FALSE]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        lag = lags[li], term = NA_character_,
        mean_beta_raw = NA_real_,
        mean_beta_corrected = NA_real_,
        mean_beta_multilag = NA_real_,
        mean_beta_unweighted = NA_real_,
        mean_beta_weighted = NA_real_,
        mean_beta_gls = NA_real_,
        mean_r_raw = NA_real_,
        mean_r_demixed = NA_real_,
        mean_tau2_bio = fmean(sub_df$tau2_bio),
        mean_gain = fmean(sub_df$mean_gain),
        mean_se_trial = fmean(sub_df$mean_se_trial),
        mean_error_ac1 = fmean(sub_df$mean_error_ac1),
        n_voxels = sum(is.finite(sub_df$tau2_bio)),
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
    cov_full <- if (n_terms > 0L) as.data.frame(do.call(rbind, cov_list)) else NULL

    panel_list <- vector("list", n_lags)
    for (li in seq_len(n_lags)) {
      base_meta <- data.table::data.table(
        run = rep(run_vec, times = n_units),
        trial = rep(trial_vec, times = n_units),
        event_index = rep(evt_idx_vec, times = n_units),
        onset = rep(onset_vec, times = n_units),
        evt_time = lags[li],
        voxel = rep(unit_ids, each = total_trials),
        decon_raw = as.numeric(raw_epochs[[li]]),
        decon_demixed = as.numeric(demixed_epochs[[li]]),
        decon_eb = as.numeric(eb_y_epochs[[li]]),
        se_raw = as.numeric(se_raw_epochs[[li]]),
        se_demixed = as.numeric(se_demixed_epochs[[li]]),
        reliability_trial = as.numeric(rel_epochs[[li]]),
        precision_weight = as.numeric(pw_epochs[[li]]),
        precision_weight_op = as.numeric(pw_op_epochs[[li]]),
        precision_weight_raw = as.numeric(pw_raw_epochs[[li]]),
        error_ac1 = as.numeric(error_ac1_epochs[[li]]),
        tau2_bio = as.numeric(tau2_bio_epochs[[li]])
      )
      if (isTRUE(aggregate_by_parcel)) {
        base_meta[["parcel"]] <- base_meta[["voxel"]]
      }
      if (!is.null(cov_full)) {
        for (cn in colnames(cov_full)) {
          base_meta[[cn]] <- rep(cov_full[[cn]], times = n_units)
        }
      }
      if (isTRUE(whiten_panel)) {
        base_meta[["decon_demixed_whitened"]] <- as.numeric(whitened_y_epochs[[li]])
        base_meta[["se_whitened"]] <- as.numeric(se_whitened_epochs[[li]])
        base_meta[["decon_eb_whitened"]] <- as.numeric(eb_whitened_y_epochs[[li]])
        if (n_terms > 0L) {
          for (ti in seq_len(n_terms)) {
            cn_w <- paste0(term_names[ti], "_whitened")
            base_meta[[cn_w]] <- as.numeric(whitened_x_epochs[[li]][[ti]])
          }
        }
      }
      panel_list[[li]] <- base_meta
    }
    panel_dt <- data.table::rbindlist(panel_list)
    attr(panel_dt, "lag_error_cov") <- lag_error_cov
    attr(panel_dt, "lag_error_cor") <- lag_error_cor
    attr(panel_dt, "lag_resolution") <- lag_resolution
  }

  list(
    coefficients = coef_df,
    summary = do.call(rbind, summary_rows),
    epochs = if (return_epochs) {
      list(
        raw = raw_epochs,
        demixed = demixed_epochs,
        eb = eb_y_epochs,
        se_raw = se_raw_epochs,
        se_demixed = se_demixed_epochs,
        reliability_trial = rel_epochs,
        precision_weight = pw_epochs,
        precision_weight_op = pw_op_epochs,
        error_ac1 = error_ac1_epochs,
        tau2_bio = tau2_bio_epochs,
        whitened = whitened_y_epochs,
        se_whitened = se_whitened_epochs,
        eb_whitened = eb_whitened_y_epochs
      )
    } else NULL,
    panel = panel_dt,
    parcel_weights = parcel_weights_list,
    lag_error_cov = lag_error_cov,
    lag_error_cor = lag_error_cor,
    lag_resolution = lag_resolution,
    settings = list(
      target_condition = target_condition,
      lags = lags,
      readout_kernel = readout_kernel,
      readout_width = readout_width,
      demix_other_events = demix_other_events,
      demix_target_trials = demix_target_trials,
      demix_within_trial = demix_within_trial,
      demix_ridge = demix_ridge,
      within_trial_ridge = within_trial_ridge,
      gain_floor = gain_floor,
      restore_scale = restore_scale,
      source_model = source_model,
      duration_scale = duration_scale,
      weighting = weighting,
      tau2_mode = eff_tau2_mode,
      aggregate_by_parcel = aggregate_by_parcel,
      whiten_panel = whiten_panel
    )
  )
}

# Smallest finite real generalized eigenvalue solving det(M - lambda * Psi) = 0
# over the error-prone subspace (projecting out error-free zero columns of Psi via Schur complement).
medusa_smallest_det_root <- function(M, Psi) {
  if (requireNamespace("geigen", quietly = TRUE)) {
    roots <- tryCatch(
      geigen::geigen(M, Psi, symmetric = FALSE, only.values = TRUE)$values,
      error = function(e) NULL
    )
    if (!is.null(roots)) {
      if (is.complex(roots)) {
        sc <- suppressWarnings(max(abs(roots), na.rm = TRUE))
        if (!is.finite(sc)) sc <- 1
        roots <- ifelse(abs(Im(roots)) <= sqrt(.Machine$double.eps) * max(1, sc), Re(roots), NA_real_)
      }
      roots <- roots[is.finite(roots)]
      if (length(roots)) return(min(roots))
    }
  }
  d_psi <- diag(Psi)
  E_idx <- which(is.finite(d_psi) & d_psi > 1e-14)
  if (!length(E_idx)) return(NA_real_)
  F_idx <- setdiff(seq_len(nrow(M)), E_idx)
  M_EF <- if (length(F_idx)) {
    M_FF <- M[F_idx, F_idx, drop = FALSE]
    M[E_idx, E_idx, drop = FALSE] -
      M[E_idx, F_idx, drop = FALSE] %*% solve(M_FF + 1e-12 * diag(length(F_idx)), M[F_idx, E_idx, drop = FALSE])
  } else {
    M[E_idx, E_idx, drop = FALSE]
  }
  Psi_EE <- 0.5 * (Psi[E_idx, E_idx, drop = FALSE] + t(Psi[E_idx, E_idx, drop = FALSE]))
  eig_P <- eigen(Psi_EE, symmetric = TRUE)
  pos <- eig_P$values > 1e-12 * max(1, max(eig_P$values))
  if (!any(pos)) return(NA_real_)
  U_pos <- eig_P$vectors[, pos, drop = FALSE]
  d_inv_half <- 1 / sqrt(eig_P$values[pos])
  H <- sweep(sweep(crossprod(U_pos, M_EF %*% U_pos), 1L, d_inv_half, "*"), 2L, d_inv_half, "*")
  H <- 0.5 * (H + t(H))
  vals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  vals <- vals[is.finite(vals)]
  if (length(vals)) min(vals) else NA_real_
}

medusa_mmom_correction <- function(lambda, nobs, alpha) {
  if (is.finite(lambda) && lambda <= 1 + 1 / nobs) {
    lambda - (1 + alpha) / nobs
  } else {
    1 - alpha / nobs
  }
}

# Built-in fallback implementation of Vig & Hallquist (2026) M-MOM / M-MOMS
# matching repair::fit_mmom() when the external `repair` package is not installed.
medusa_fit_mmom_fallback <- function(
  data_list,
  preliminary_moment_control = list(correction = "step_up", alpha = NULL),
  final_moment_control = list(correction = "step_up", alpha = NULL)
) {
  y <- as.numeric(data_list$outcome)
  X_pred <- as.matrix(data_list$predictors)
  storage.mode(X_pred) <- "double"
  nobs <- length(y)
  q <- ncol(X_pred)
  p <- q + 1L
  if (is.null(colnames(X_pred))) colnames(X_pred) <- paste0("x", seq_len(q))
  design <- cbind(`(Intercept)` = 1.0, X_pred)
  pred_idx <- 2L:p
  err_x <- data_list$err_x
  err_y <- if (is.null(data_list$err_y)) rep(0.0, nobs) else as.numeric(data_list$err_y)
  err_xy <- if (is.null(data_list$err_xy)) matrix(0.0, nobs, q) else as.matrix(data_list$err_xy)

  alpha_min <- q + 2.0
  norm_ctrl <- function(ctrl) {
    corr <- if (is.null(ctrl$correction)) "minimum" else ctrl$correction
    alp <- if (corr == "custom" && !is.null(ctrl$alpha)) as.numeric(ctrl$alpha) else alpha_min
    list(correction = corr, alpha = alp)
  }
  c_pre <- norm_ctrl(preliminary_moment_control)
  c_fin <- norm_ctrl(final_moment_control)

  sum_err_x <- function(w) {
    out <- matrix(0.0, q, q)
    for (i in seq_len(nobs)) {
      out <- out + w[i] * err_x[, , i]
    }
    0.5 * (out + t(out))
  }
  build_err_moment <- function(ex_sum, ey_sum, exy_sum_design) {
    em <- matrix(0.0, p + 1L, p + 1L)
    em[1L, 1L] <- ey_sum
    em[1L, 2L:(p + 1L)] <- exy_sum_design
    em[2L:(p + 1L), 1L] <- exy_sum_design
    em[2L:(p + 1L), 2L:(p + 1L)] <- ex_sum
    em
  }

  ex_sum_1 <- sum_err_x(rep(1.0, nobs))
  err_cov_design_1 <- matrix(0.0, p, p)
  err_cov_design_1[pred_idx, pred_idx] <- ex_sum_1
  exy_sum_1 <- c(0.0, colSums(err_xy))
  err_mom_1 <- build_err_moment(err_cov_design_1, sum(err_y), exy_sum_1)
  obs_mom_1 <- crossprod(cbind(y, design))
  lam_1 <- medusa_smallest_det_root(obs_mom_1, err_mom_1)
  des_mom_1 <- crossprod(design)
  out_mom_1 <- as.numeric(crossprod(design, y))

  eval_fixed <- function(alpha_1, corr_1, alpha_2, corr_2) {
    c1 <- if (identical(corr_1, "none") || max(abs(err_mom_1)) <= .Machine$double.eps) {
      1.0
    } else {
      max(0.0, medusa_mmom_correction(lam_1, nobs, alpha_1))
    }
    S1_x <- des_mom_1 - c1 * err_cov_design_1
    s1_xy <- out_mom_1 - c1 * exy_sum_1
    coef0 <- tryCatch(as.numeric(qr.solve(S1_x, s1_xy)), error = function(e) NULL)
    if (is.null(coef0) || any(!is.finite(coef0))) {
      return(list(converged = FALSE, S1_x = S1_x, c1 = c1, lam_1 = lam_1))
    }
    names(coef0) <- colnames(design)
    g0_pred <- coef0[pred_idx]
    err_resid_var <- numeric(nobs)
    for (i in seq_len(nobs)) {
      err_resid_var[i] <- err_y[i] - 2.0 * sum(g0_pred * err_xy[i, ]) +
        as.numeric(crossprod(g0_pred, err_x[, , i] %*% g0_pred))
    }
    err_resid_var <- pmax(0.0, err_resid_var)
    resid1 <- y - as.numeric(design %*% coef0)
    sigma2_eps <- if (is.finite(lam_1) && lam_1 < 1.0) {
      0.0
    } else {
      max(0.0, sum(resid1^2) / max(1L, nobs - p) - mean(err_resid_var))
    }
    case_var <- sigma2_eps + err_resid_var
    if (any(!is.finite(case_var)) || any(case_var <= sqrt(.Machine$double.eps))) {
      return(list(converged = FALSE, S1_x = S1_x, c1 = c1, lam_1 = lam_1))
    }
    w_case <- 1.0 / case_var
    w_des <- design * sqrt(w_case)
    w_out <- y * sqrt(w_case)
    w_mom <- crossprod(cbind(w_out, w_des))
    ex_sum_2 <- sum_err_x(w_case)
    err_cov_design_2 <- matrix(0.0, p, p)
    err_cov_design_2[pred_idx, pred_idx] <- ex_sum_2
    exy_sum_2 <- c(0.0, colSums(err_xy * w_case))
    err_mom_2 <- build_err_moment(err_cov_design_2, sum(w_case * err_y), exy_sum_2)
    lam_2 <- medusa_smallest_det_root(w_mom, err_mom_2)
    c2 <- if (identical(corr_2, "none") || max(abs(err_mom_2)) <= .Machine$double.eps) {
      1.0
    } else {
      max(0.0, medusa_mmom_correction(lam_2, nobs, alpha_2))
    }
    corr_mom_2 <- w_mom - c2 * err_mom_2
    S2_x <- corr_mom_2[2L:(p + 1L), 2L:(p + 1L), drop = FALSE]
    s2_xy <- corr_mom_2[2L:(p + 1L), 1L]
    coef2 <- tryCatch(as.numeric(qr.solve(S2_x, s2_xy)), error = function(e) NULL)
    if (is.null(coef2) || any(!is.finite(coef2))) {
      return(list(converged = FALSE, S1_x = S1_x, S2_x = S2_x, c1 = c1, c2 = c2, lam_1 = lam_1, lam_2 = lam_2))
    }
    names(coef2) <- colnames(design)
    bread_inv <- tryCatch(qr.solve(S2_x / nobs, diag(p)), error = function(e) NULL)
    if (is.null(bread_inv) || any(!is.finite(bread_inv))) {
      return(list(converged = FALSE, S1_x = S1_x, S2_x = S2_x, c1 = c1, c2 = c2, lam_1 = lam_1, lam_2 = lam_2))
    }
    tilde_cp <- matrix(0.0, p, p)
    for (i in seq_len(nobs)) {
      om_i <- c(0.0, err_xy[i, ] - as.numeric(err_x[, , i] %*% g0_pred))
      tilde_cp <- tilde_cp + (w_case[i]^2) * tcrossprod(om_i)
    }
    meat <- crossprod(w_des) + tilde_cp
    vcov_mat <- 0.5 * ((bread_inv %*% meat %*% bread_inv) / (nobs^2) +
      t((bread_inv %*% meat %*% bread_inv) / (nobs^2)))
    dimnames(vcov_mat) <- list(colnames(design), colnames(design))
    if (any(!is.finite(diag(vcov_mat))) || any(diag(vcov_mat) < 0)) {
      return(list(converged = FALSE, S1_x = S1_x, S2_x = S2_x, c1 = c1, c2 = c2, lam_1 = lam_1, lam_2 = lam_2))
    }
    se_vec <- sqrt(diag(vcov_mat))
    z_vec <- coef2 / se_vec
    p_vec <- 2.0 * stats::pnorm(abs(z_vec), lower.tail = FALSE)
    coef_df <- data.frame(
      term = colnames(design),
      estimate = coef2,
      se = se_vec,
      zstat = z_vec,
      pval = p_vec,
      ci_low = coef2 - 1.959964 * se_vec,
      ci_high = coef2 + 1.959964 * se_vec,
      stringsAsFactors = FALSE
    )
    list(
      converged = TRUE,
      coefficients = coef_df,
      vcov = vcov_mat,
      preliminary_coef = coef0,
      corrected_design_moment = S2_x,
      analysis_weights = w_case,
      S1_x = S1_x,
      S2_x = S2_x,
      w_des_mom = crossprod(w_des),
      c1 = c1,
      c2 = c2,
      lam_1 = lam_1,
      lam_2 = lam_2,
      alpha_prelim = alpha_1,
      alpha_final = alpha_2,
      resid_var = sigma2_eps
    )
  }

  check_guard <- function(fit_cand, stage) {
    if (!isTRUE(fit_cand$converged)) return(FALSE)
    if (stage == "preliminary") {
      S_corr <- fit_cand$S1_x[pred_idx, pred_idx, drop = FALSE] / nobs
      S_obs <- des_mom_1[pred_idx, pred_idx, drop = FALSE] / nobs
    } else {
      S_corr <- fit_cand$S2_x[pred_idx, pred_idx, drop = FALSE]
      S_obs <- fit_cand$w_des_mom[pred_idx, pred_idx, drop = FALSE]
    }
    eig_c <- eigen(0.5 * (S_corr + t(S_corr)), symmetric = TRUE, only.values = TRUE)$values
    eig_o <- eigen(0.5 * (S_obs + t(S_obs)), symmetric = TRUE, only.values = TRUE)$values
    min_c <- min(eig_c)
    max_c <- max(eig_c)
    max_o <- max(eig_o)
    is.finite(min_c) && is.finite(max_c) && is.finite(max_o) &&
      (min_c >= sqrt(.Machine$double.eps)) &&
      (min_c >= 0.05 * max_o) &&
      ((max_c / max(min_c, .Machine$double.eps)) <= 1e5)
  }

  alpha_1 <- c_pre$alpha
  alpha_2 <- c_fin$alpha
  fit_cur <- eval_fixed(
    alpha_1, if (c_pre$correction == "step_up") "custom" else c_pre$correction,
    alpha_2, if (c_fin$correction == "step_up") "custom" else c_fin$correction
  )

  if (c_pre$correction == "step_up" && !check_guard(fit_cur, "preliminary")) {
    upper_1 <- if (is.finite(lam_1)) {
      if (lam_1 <= 1 + 1 / nobs) lam_1 * nobs - 1 else nobs
    } else alpha_min
    if (upper_1 > alpha_min) {
      grid_1 <- seq(alpha_min, upper_1, length.out = 9L)
      pass_hi <- NA_real_
      fail_lo <- alpha_min
      fit_hi <- NULL
      for (a_cand in grid_1[-1L]) {
        f_c <- eval_fixed(a_cand, "custom", alpha_2, if (c_fin$correction == "step_up") "custom" else c_fin$correction)
        if (check_guard(f_c, "preliminary")) {
          pass_hi <- a_cand
          fit_hi <- f_c
          break
        } else {
          fail_lo <- a_cand
        }
      }
      if (is.finite(pass_hi)) {
        lo <- fail_lo; hi <- pass_hi
        for (iter in seq_len(12L)) {
          if (hi - lo <= 0.005) break
          mid <- 0.5 * (lo + hi)
          f_m <- eval_fixed(mid, "custom", alpha_2, if (c_fin$correction == "step_up") "custom" else c_fin$correction)
          if (check_guard(f_m, "preliminary")) {
            hi <- mid; fit_hi <- f_m
          } else {
            lo <- mid
          }
        }
        alpha_1 <- hi
        fit_cur <- fit_hi
      }
    }
  }

  if (c_fin$correction == "step_up" && !check_guard(fit_cur, "final")) {
    lam_2_cur <- fit_cur$lam_2
    upper_2 <- if (!is.null(lam_2_cur) && is.finite(lam_2_cur)) {
      if (lam_2_cur <= 1 + 1 / nobs) lam_2_cur * nobs - 1 else nobs
    } else alpha_min
    if (upper_2 > alpha_min) {
      grid_2 <- seq(alpha_min, upper_2, length.out = 9L)
      pass_hi <- NA_real_
      fail_lo <- alpha_min
      fit_hi <- NULL
      for (a_cand in grid_2[-1L]) {
        f_c <- eval_fixed(alpha_1, if (c_pre$correction == "step_up") "custom" else c_pre$correction, a_cand, "custom")
        if (check_guard(f_c, "final")) {
          pass_hi <- a_cand
          fit_hi <- f_c
          break
        } else {
          fail_lo <- a_cand
        }
      }
      if (is.finite(pass_hi)) {
        lo <- fail_lo; hi <- pass_hi
        for (iter in seq_len(12L)) {
          if (hi - lo <= 0.005) break
          mid <- 0.5 * (lo + hi)
          f_m <- eval_fixed(alpha_1, if (c_pre$correction == "step_up") "custom" else c_pre$correction, mid, "custom")
          if (check_guard(f_m, "final")) {
            hi <- mid; fit_hi <- f_m
          } else {
            lo <- mid
          }
        }
        alpha_2 <- hi
        fit_cur <- fit_hi
      }
    }
  }

  fit_cur$diagnostics <- list(
    nobs = nobs,
    npredictors = q,
    lambda_prelim = fit_cur$lam_1,
    lambda_final = fit_cur$lam_2,
    alpha_prelim = alpha_1,
    alpha_final = alpha_2,
    correction_prelim = fit_cur$c1,
    correction_final = fit_cur$c2,
    resid_var = fit_cur$resid_var
  )
  fit_cur
}

#' Fit a REPAIR modified method-of-moments (M-MOM / M-MOMS) measurement-error model
#'
#' Estimates a linear regression model when predictors, the outcome, or both are
#' likelihood-only predictions measured with heteroscedastic sampling error
#' (Vig & Hallquist, 2026). Delegates directly to \code{repair::fit_mmom()} when
#' the \pkg{repair} package is installed and uses an identical built-in fallback
#' otherwise.
#'
#' @param outcome Numeric outcome vector of length \eqn{N}.
#' @param predictors Numeric predictor matrix (\eqn{N \times q}, excluding the
#'   intercept).
#' @param err_x Numeric 3D array (\eqn{q \times q \times N}) of per-observation
#'   predictor sampling-error covariance matrices \eqn{\hat{\boldsymbol{\Psi}}_{x,j}}.
#' @param err_y Non-negative numeric vector of length \eqn{N} giving per-observation
#'   outcome sampling-error variances \eqn{\hat{\psi}_{y,j}}.
#' @param err_xy Optional numeric matrix (\eqn{N \times q}) of per-observation
#'   outcome-predictor sampling-error covariances \eqn{\hat{\boldsymbol{\psi}}_{xy,j}^\top}.
#' @param regularization Regularization policy: \code{"step_up"} (default;
#'   adaptive M-MOMS step-up search), \code{"minimum"} (fixed minimum
#'   \eqn{\alpha = q + 2}), \code{"none"}, or \code{"custom"}.
#' @param alpha Optional custom regularization parameter when
#'   \code{regularization = "custom"}.
#' @param backend Computation backend passed to \code{repair::fit_mmom()} when
#'   available (\code{"auto"}, \code{"r"}, or \code{"rcpp"}).
#' @return A list with \code{converged}, \code{coefficients} (data.frame with
#'   \code{term}, \code{estimate}, \code{se}, \code{zstat}, \code{pval},
#'   \code{ci_low}, \code{ci_high}), \code{vcov}, \code{preliminary_coef},
#'   \code{corrected_design_moment}, \code{analysis_weights},
#'   \code{diagnostics}, and \code{repair_fit} (the raw \pkg{repair} object when
#'   \pkg{repair} is installed).
#' @export
medusa_fit_mmom <- function(
  outcome,
  predictors,
  err_x,
  err_y = NULL,
  err_xy = NULL,
  regularization = c("step_up", "minimum", "none", "custom"),
  alpha = NULL,
  backend = c("auto", "r", "rcpp")
) {
  regularization <- match.arg(regularization)
  backend <- match.arg(backend)
  predictors <- as.matrix(predictors)
  storage.mode(predictors) <- "double"
  nobs <- length(outcome)
  q <- ncol(predictors)
  if (is.null(colnames(predictors))) colnames(predictors) <- paste0("x", seq_len(q))
  if (is.null(err_y)) err_y <- rep(0.0, nobs)
  if (is.null(err_xy)) err_xy <- matrix(0.0, nobs, q)
  ctrl <- list(correction = regularization, alpha = alpha)
  data_list <- list(
    outcome = as.numeric(outcome),
    predictors = predictors,
    err_x = err_x,
    err_y = as.numeric(err_y),
    err_xy = as.matrix(err_xy)
  )

  if (requireNamespace("repair", quietly = TRUE)) {
    rfit <- tryCatch(
      repair::fit_mmom(
        data_list = data_list,
        preliminary_moment_control = ctrl,
        final_moment_control = ctrl,
        backend = backend
      ),
      error = function(e) NULL
    )
    if (!is.null(rfit) && isTRUE(rfit@converged)) {
      diag_list <- rfit@diagnostics
      return(list(
        converged = TRUE,
        coefficients = as.data.frame(rfit@coefficients),
        vcov = as.matrix(rfit@vcov),
        preliminary_coef = rfit@preliminary_coef,
        corrected_design_moment = as.matrix(rfit@corrected_design_moment),
        analysis_weights = as.numeric(rfit@analysis_weights),
        diagnostics = list(
          nobs = nobs,
          npredictors = q,
          lambda_prelim = diag_list$lambda_prelim,
          lambda_final = diag_list$lambda_final,
          alpha_prelim = if (!is.null(diag_list$alpha_prelim)) diag_list$alpha_prelim else q + 2,
          alpha_final = if (!is.null(diag_list$alpha_final)) diag_list$alpha_final else q + 2
        ),
        repair_fit = rfit
      ))
    }
  }

  fb <- medusa_fit_mmom_fallback(
    data_list = data_list,
    preliminary_moment_control = ctrl,
    final_moment_control = ctrl
  )
  fb$repair_fit <- NULL
  fb
}

#' Measurement-error-corrected inter-region connectivity and PPI via REPAIR
#'
#' Estimates trial-by-trial seed-to-target neural connectivity, directed lagged
#' connectivity (controlling for baseline target activity at the seed lag), and
#' task-modulated psychophysiological interactions (PPI) from a
#' \code{\link{decon_epoch_regression}} trial panel using the REPAIR modified
#' method-of-moments (\code{M-MOMS} / \code{M-MOM}) framework of Vig & Hallquist
#' (2026). Propagates per-trial operator standard errors (\code{se_demixed}) into
#' exact heteroscedastic predictor and outcome error covariances,including the
#' rank-1 interaction error block for PPI (\code{seed * x}) and dual predictor
#' error variances for Granger-style auto-lag controls (\code{targ_t0}).
#'
#' @param panel_dt A \code{data.frame} or \code{data.table} returned in
#'   \code{$panel} by \code{\link{decon_epoch_regression}}.
#' @param seed_region Character scalar naming the seed region/parcel/voxel.
#' @param target_region Character scalar naming the target region/parcel/voxel.
#' @param seed_lag Numeric lag (in seconds) for the seed readout (default \code{0}).
#' @param target_lag Numeric lag (in seconds) for the target readout (default \code{0}).
#' @param covariates Optional character vector of error-free trial covariate
#'   columns in \code{panel_dt}. Defaults to \code{"x"} when present in
#'   \code{panel_dt}, or \code{character(0)} otherwise.
#' @param test_ppi Logical; if \code{TRUE}, includes the interaction between
#'   the seed readout and the first covariate in \code{covariates}
#'   (\code{seed_c:x_c}) with exact per-trial interaction error covariance
#'   \eqn{\tilde{\psi}_{S,j} [1, x_{c,j}; x_{c,j}, x_{c,j}^2]}. Default: \code{FALSE}.
#' @param control_target_seed_lag Logical; if \code{TRUE} and
#'   \code{abs(target_lag - seed_lag) > 1e-6}, includes the target region's
#'   readout at \code{seed_lag} (\code{targ_t0_val}) as an error-prone control
#'   covariate with measurement-error variance \code{se_demixed^2}. Default:
#'   \code{FALSE}.
#' @param orthogonalize_seed_x Logical; if \code{TRUE} (default), residualizes
#'   the run-centered seed readout on the error-free covariates within each
#'   subject-run and scales the seed error variance by \eqn{1 - h_{x,r,j}}
#'   prior to forming PPI interactions, preventing quadratic task drive
#'   (\eqn{x_c^2}) from confounding physiological-by-psychological interactions.
#' @param method Estimation method: \code{"mmoms"} (default; adaptive step-up
#'   REPAIR \code{M-MOMS}), \code{"mmom"} (fixed minimum-regularization REPAIR
#'   \code{M-MOM}), \code{"eb_predictor"} (Skrondal & Laake 2001 role-specific
#'   Empirical-Bayes shrunk predictor on likelihood-only outcome),
#'   \code{"naive_gls"} (precision-weighted WLS without predictor error
#'   subtraction), or \code{"naive_unweighted"} (unweighted OLS).
#' @param use_whitened Logical; if \code{TRUE}, uses \code{decon_demixed_whitened}
#'   and \code{se_whitened} when present. Defaults to \code{FALSE} for
#'   \code{"mmoms"} and \code{"mmom"} (so elementwise PPI products are formed on
#'   unmixed trial states before error correction) and \code{TRUE} for
#'   \code{"naive_gls"}.
#' @param region_col Optional column name identifying regions/parcels/voxels in
#'   \code{panel_dt} (auto-detected from \code{"region"}, \code{"parcel"},
#'   \code{"voxel"}).
#' @param id_col Column name identifying subjects in \code{panel_dt} (default
#'   \code{"id"}).
#' @param return_augmented_panel Logical; if \code{TRUE}, returns the merged
#'   trial panel augmented with moment-preserving REPAIR predictor columns
#'   (\code{seed_repair}, \code{ppi_repair}, \code{repair_weight}) that reproduce
#'   subject-level \code{M-MOMS} slopes in weighted least squares / \code{lmer()}.
#' @return A list of class \code{"medusa_connectivity_repair"} containing
#'   \code{group} (group-level random-effects meta-analytic inference table),
#'   \code{subjects} (per-subject estimates, SEs, and REPAIR diagnostics),
#'   optional \code{panel_repair}, and \code{settings}.
#' @export
decon_connectivity_repair <- function(
  panel_dt,
  seed_region,
  target_region,
  seed_lag = 0,
  target_lag = 0,
  covariates = NULL,
  test_ppi = FALSE,
  control_target_seed_lag = FALSE,
  orthogonalize_seed_x = TRUE,
  method = c("mmoms", "mmom", "eb_predictor", "naive_gls", "naive_unweighted"),
  use_whitened = NULL,
  region_col = NULL,
  id_col = "id",
  return_augmented_panel = FALSE
) {
  method <- match.arg(method)
  dt <- data.table::copy(data.table::as.data.table(panel_dt))
  if (is.null(region_col)) {
    for (cand in c("region", "parcel", "voxel")) {
      if (cand %in% names(dt)) {
        region_col <- cand
        break
      }
    }
  }
  if (is.null(region_col) || !(region_col %in% names(dt))) {
    stop("Could not find region/parcel/voxel column in panel_dt")
  }
  if (!(id_col %in% names(dt))) {
    dt[, (id_col) := "sub001"]
  }
  if (!("run" %in% names(dt))) {
    dt[, run := "run1"]
  }
  if (is.null(covariates)) {
    covariates <- if ("x" %in% names(dt)) "x" else character(0L)
  } else {
    covariates <- intersect(covariates, names(dt))
  }
  if (isTRUE(test_ppi) && !length(covariates)) {
    stop("test_ppi = TRUE requires at least one covariate column in panel_dt")
  }
  if (is.null(use_whitened)) {
    use_whitened <- identical(method, "naive_gls")
  }

  y_col <- if (isTRUE(use_whitened) && "decon_demixed_whitened" %in% names(dt)) {
    "decon_demixed_whitened"
  } else if ("decon_demixed" %in% names(dt)) {
    "decon_demixed"
  } else if ("beta_lss" %in% names(dt)) {
    "beta_lss"
  } else {
    "decon_raw"
  }
  se_col <- if (isTRUE(use_whitened) && "se_whitened" %in% names(dt)) {
    "se_whitened"
  } else if ("se_demixed" %in% names(dt)) {
    "se_demixed"
  } else {
    NULL
  }
  eb_col <- if (isTRUE(use_whitened) && "decon_eb_whitened" %in% names(dt)) {
    "decon_eb_whitened"
  } else if ("decon_eb" %in% names(dt)) {
    "decon_eb"
  } else {
    NULL
  }

  seed_sub <- dt[get(region_col) == seed_region & abs(evt_time - seed_lag) < 1e-6]
  targ_sub <- dt[get(region_col) == target_region & abs(evt_time - target_lag) < 1e-6]
  if (!nrow(seed_sub)) stop("No rows match seed_region = '", seed_region, "' at seed_lag = ", seed_lag)
  if (!nrow(targ_sub)) stop("No rows match target_region = '", target_region, "' at target_lag = ", target_lag)

  seed_keep <- data.table::data.table(
    id = seed_sub[[id_col]],
    run = seed_sub[["run"]],
    trial = seed_sub[["trial"]],
    seed_val = as.numeric(seed_sub[[y_col]]),
    se_seed = if (!is.null(se_col)) as.numeric(seed_sub[[se_col]]) else rep(0.0, nrow(seed_sub)),
    seed_eb = if (!is.null(eb_col)) as.numeric(seed_sub[[eb_col]]) else NA_real_,
    tau2_seed = if ("tau2_bio" %in% names(seed_sub)) as.numeric(seed_sub[["tau2_bio"]]) else NA_real_
  )
  for (cv in covariates) {
    cv_use <- if (isTRUE(use_whitened) && paste0(cv, "_whitened") %in% names(seed_sub)) {
      paste0(cv, "_whitened")
    } else {
      cv
    }
    seed_keep[[cv]] <- as.numeric(seed_sub[[cv_use]])
  }

  targ_keep <- data.table::data.table(
    id = targ_sub[[id_col]],
    run = targ_sub[["run"]],
    trial = targ_sub[["trial"]],
    targ_val = as.numeric(targ_sub[[y_col]]),
    se_targ = if (!is.null(se_col)) as.numeric(targ_sub[[se_col]]) else rep(0.0, nrow(targ_sub)),
    w_targ = if ("precision_weight" %in% names(targ_sub)) as.numeric(targ_sub[["precision_weight"]]) else rep(1.0, nrow(targ_sub))
  )

  m_dt <- merge(targ_keep, seed_keep, by = c("id", "run", "trial"), sort = FALSE)
  has_t0 <- isTRUE(control_target_seed_lag) && abs(target_lag - seed_lag) > 1e-6
  if (has_t0) {
    t0_sub <- dt[get(region_col) == target_region & abs(evt_time - seed_lag) < 1e-6]
    if (nrow(t0_sub) > 0L) {
      t0_keep <- data.table::data.table(
        id = t0_sub[[id_col]],
        run = t0_sub[["run"]],
        trial = t0_sub[["trial"]],
        targ_t0_val = as.numeric(t0_sub[[y_col]]),
        se_targ_t0 = if (!is.null(se_col)) as.numeric(t0_sub[[se_col]]) else rep(0.0, nrow(t0_sub)),
        targ_t0_eb = if (!is.null(eb_col)) as.numeric(t0_sub[[eb_col]]) else NA_real_,
        tau2_targ_t0 = if ("tau2_bio" %in% names(t0_sub)) as.numeric(t0_sub[["tau2_bio"]]) else NA_real_
      )
      m_dt <- merge(m_dt, t0_keep, by = c("id", "run", "trial"), sort = FALSE)
    } else {
      has_t0 <- FALSE
    }
  }

  # Center target, seed, and covariates within (id, run)
  m_dt[, targ_c := targ_val - mean(targ_val, na.rm = TRUE), by = .(id, run)]
  m_dt[, seed_c_raw := seed_val - mean(seed_val, na.rm = TRUE), by = .(id, run)]
  if (has_t0) {
    m_dt[, targ_t0_c := targ_t0_val - mean(targ_t0_val, na.rm = TRUE), by = .(id, run)]
  }
  for (cv in covariates) {
    m_dt[, (paste0(cv, "_c")) := get(cv) - mean(get(cv), na.rm = TRUE), by = .(id, run)]
  }
  cov_c_names <- if (length(covariates)) paste0(covariates, "_c") else character(0L)

  # Orthogonalize seed_c with respect to error-free covariates within (id, run) when requested
  if (isTRUE(orthogonalize_seed_x) && length(cov_c_names) > 0L) {
    m_dt[, c("seed_c", "h_x_run") := {
      X_r <- as.matrix(.SD)
      ys <- seed_c_raw
      ok_r <- is.finite(ys) & (rowSums(!is.finite(X_r)) == 0L)
      if (sum(ok_r) > ncol(X_r) + 2L && any(apply(X_r[ok_r, , drop = FALSE], 2L, stats::sd) > 1e-10)) {
        qr_x <- qr(X_r[ok_r, , drop = FALSE])
        Q_x <- qr.Q(qr_x)
        res_s <- ys
        res_s[ok_r] <- ys[ok_r] - as.numeric(Q_x %*% crossprod(Q_x, ys[ok_r]))
        h_vec <- rep(0.0, length(ys))
        h_vec[ok_r] <- pmin(0.85, rowSums(Q_x^2))
        list(res_s, h_vec)
      } else {
        list(ys, rep(0.0, length(ys)))
      }
    }, by = .(id, run), .SDcols = cov_c_names]
  } else {
    m_dt[, seed_c := seed_c_raw]
    m_dt[, h_x_run := 0.0]
  }

  # Compute EB reliability factor per trial for role-specific EB predictor mode
  m_dt[, rel_seed_trial := {
    s2 <- (1.0 - h_x_run) * (se_seed^2)
    ok_s <- is.finite(s2) & (s2 > 0)
    if (any(ok_s)) {
      v_obs <- stats::var(seed_c, na.rm = TRUE)
      t2_est <- max(0.05 * stats::median(s2[ok_s]), v_obs - mean(s2[ok_s]))
      t2_est / (t2_est + s2)
    } else {
      rep(1.0, .N)
    }
  }, by = id]
  m_dt[, seed_eb_c := rel_seed_trial * seed_c]
  if (has_t0) {
    m_dt[, rel_t0_trial := {
      s2 <- se_targ_t0^2
      ok_s <- is.finite(s2) & (s2 > 0)
      if (any(ok_s)) {
        v_obs <- stats::var(targ_t0_c, na.rm = TRUE)
        t2_est <- max(0.05 * stats::median(s2[ok_s]), v_obs - mean(s2[ok_s]))
        t2_est / (t2_est + s2)
      } else {
        rep(1.0, .N)
      }
    }, by = id]
    m_dt[, targ_t0_eb_c := rel_t0_trial * targ_t0_c]
  }

  subj_ids <- unique(m_dt$id)
  subj_rows <- vector("list", length(subj_ids))
  aug_panels <- if (isTRUE(return_augmented_panel)) vector("list", length(subj_ids)) else NULL

  for (si in seq_along(subj_ids)) {
    sid <- subj_ids[si]
    sub_dt <- m_dt[id == sid]
    n_i <- nrow(sub_dt)

    use_eb <- identical(method, "eb_predictor")
    s_vec <- if (use_eb) sub_dt$seed_eb_c else sub_dt$seed_c
    psi_s <- if (method %in% c("mmoms", "mmom")) {
      pmax(0.0, (1.0 - sub_dt$h_x_run) * (sub_dt$se_seed^2))
    } else {
      rep(0.0, n_i)
    }
    psi_y <- if (method %in% c("mmoms", "mmom")) {
      pmax(0.0, sub_dt$se_targ^2)
    } else {
      rep(0.0, n_i)
    }

    # Build predictor matrix X_i and per-trial error covariance array err_x_i
    pred_list <- list(seed_c = s_vec)
    if (has_t0) {
      pred_list[["targ_t0_val"]] <- if (use_eb) sub_dt$targ_t0_eb_c else sub_dt$targ_t0_c
    }
    for (k in seq_along(covariates)) {
      pred_list[[ paste0(covariates[k], "_c") ]] <- sub_dt[[ cov_c_names[k] ]]
    }
    if (isTRUE(test_ppi)) {
      x1_c <- sub_dt[[ cov_c_names[1L] ]]
      pred_list[["seed_c:x_c"]] <- s_vec * x1_c
    }
    X_i <- do.call(cbind, pred_list)
    q_i <- ncol(X_i)
    err_x_i <- array(0.0, dim = c(q_i, q_i, n_i))
    err_x_i[1L, 1L, ] <- psi_s
    if (has_t0 && method %in% c("mmoms", "mmom")) {
      err_x_i[2L, 2L, ] <- pmax(0.0, sub_dt$se_targ_t0^2)
    }
    if (isTRUE(test_ppi) && method %in% c("mmoms", "mmom")) {
      x1_c <- sub_dt[[ cov_c_names[1L] ]]
      idx_ppi <- q_i
      err_x_i[1L, idx_ppi, ] <- x1_c * psi_s
      err_x_i[idx_ppi, 1L, ] <- x1_c * psi_s
      err_x_i[idx_ppi, idx_ppi, ] <- (x1_c^2) * psi_s
    }

    if (method %in% c("mmoms", "mmom")) {
      reg_pol <- if (identical(method, "mmoms")) "step_up" else "minimum"
      mfit <- medusa_fit_mmom(
        outcome = sub_dt$targ_c,
        predictors = X_i,
        err_x = err_x_i,
        err_y = psi_y,
        regularization = reg_pol
      )
      if (isTRUE(mfit$converged)) {
        cf <- mfit$coefficients
        cf <- cf[cf$term != "(Intercept)", , drop = FALSE]
        w_raw <- mfit$analysis_weights
        w_rep <- w_raw / mean(w_raw)
        S2_full <- mfit$corrected_design_moment
        S2_pred <- S2_full[-1L, -1L, drop = FALSE] -
          tcrossprod(S2_full[-1L, 1L], S2_full[1L, -1L]) / max(S2_full[1L, 1L], .Machine$double.eps)
        mu_X_w <- colSums(X_i * w_raw) / sum(w_raw)
        X_i_wcent <- sweep(X_i, 2L, mu_X_w, "-")
        X_w_obs <- crossprod(X_i_wcent, X_i_wcent * w_raw)
        R_mat <- tryCatch(qr.solve(X_w_obs, S2_pred), error = function(e) diag(q_i))
        X_rep <- X_i_wcent %*% R_mat
        colnames(X_rep) <- paste0(colnames(X_i), "_repair")
        subj_rows[[si]] <- data.table::data.table(
          id = sid,
          term = cf$term,
          estimate = cf$estimate,
          std.error = cf$se,
          statistic = cf$zstat,
          p.value = cf$pval,
          n_trials = n_i,
          lambda_prelim = mfit$diagnostics$lambda_prelim,
          lambda_final = mfit$diagnostics$lambda_final,
          alpha_prelim = mfit$diagnostics$alpha_prelim,
          alpha_final = mfit$diagnostics$alpha_final,
          converged = TRUE
        )
        if (isTRUE(return_augmented_panel)) {
          sub_copy <- data.table::copy(sub_dt)
          for (col_j in seq_len(q_i)) {
            data.table::set(sub_copy, j = colnames(X_rep)[col_j], value = X_rep[, col_j])
          }
          data.table::set(sub_copy, j = "repair_weight", value = w_rep)
          aug_panels[[si]] <- sub_copy
        }
        next
      }
    }

    # Fallback / WLS / EB predictor / OLS path
    w_use <- if (identical(method, "naive_unweighted")) {
      rep(1.0, n_i)
    } else {
      w0 <- sub_dt$w_targ
      w0[!is.finite(w0) | w0 <= 0] <- 1.0
      w0 / mean(w0)
    }
    df_lm <- as.data.frame(X_i)
    df_lm$targ_c <- sub_dt$targ_c
    lm_fit <- stats::lm(targ_c ~ ., data = df_lm, weights = w_use)
    sm_cf <- summary(lm_fit)$coefficients
    keep_rows <- setdiff(rownames(sm_cf), "(Intercept)")
    clean_terms <- gsub("^`|`$", "", keep_rows)
    subj_rows[[si]] <- data.table::data.table(
      id = sid,
      term = clean_terms,
      estimate = as.numeric(sm_cf[keep_rows, "Estimate"]),
      std.error = as.numeric(sm_cf[keep_rows, "Std. Error"]),
      statistic = as.numeric(sm_cf[keep_rows, "t value"]),
      p.value = as.numeric(sm_cf[keep_rows, "Pr(>|t|)"]),
      n_trials = n_i,
      lambda_prelim = NA_real_,
      lambda_final = NA_real_,
      alpha_prelim = NA_real_,
      alpha_final = NA_real_,
      converged = TRUE
    )
    if (isTRUE(return_augmented_panel)) {
      sub_copy <- data.table::copy(sub_dt)
      for (col_j in seq_len(q_i)) {
        data.table::set(sub_copy, j = paste0(colnames(X_i)[col_j], "_repair"), value = X_i[, col_j])
      }
      data.table::set(sub_copy, j = "repair_weight", value = w_use)
      aug_panels[[si]] <- sub_copy
    }
  }

  subj_dt <- data.table::rbindlist(subj_rows, fill = TRUE)

  # Group-level random-effects pooling with ancillary sample-size weights and HKSJ t_{S-1} calibration
  # (avoids the downward covariance bias Cov(w_s, theta_s) < 0 caused by EIV sandwich SEs depending on theta_s^2)
  terms_u <- unique(subj_dt$term)
  grp_rows <- vector("list", length(terms_u))
  for (ti in seq_along(terms_u)) {
    tm <- terms_u[ti]
    st <- subj_dt[term == tm & is.finite(estimate) & is.finite(std.error) & std.error > 0]
    S_k <- nrow(st)
    if (S_k == 0L) next
    if (S_k == 1L) {
      grp_rows[[ti]] <- data.table::data.table(
        term = tm,
        estimate = st$estimate[1L],
        std.error = st$std.error[1L],
        statistic = st$statistic[1L],
        df = Inf,
        p.value = st$p.value[1L],
        conf.low = st$estimate[1L] - 1.959964 * st$std.error[1L],
        conf.high = st$estimate[1L] + 1.959964 * st$std.error[1L],
        tau2_between = 0.0,
        mean_unweighted = st$estimate[1L],
        median_subj = st$estimate[1L],
        n_subjects = 1L
      )
      next
    }
    theta_raw <- st$estimate
    if (S_k >= 6L) {
      med_t <- stats::median(theta_raw)
      mad_t <- stats::mad(theta_raw, constant = 1.4826)
      if (is.finite(mad_t) && mad_t > 1e-8) {
        theta_s <- pmin(med_t + 3.0 * mad_t, pmax(med_t - 3.0 * mad_t, theta_raw))
      } else {
        theta_s <- theta_raw
      }
    } else {
      theta_s <- theta_raw
    }
    v_s <- st$std.error^2
    med_v <- stats::median(v_s)
    v_s_stab <- pmax(v_s, 0.15 * med_v)
    w0 <- 1.0 / v_s_stab
    mu0 <- sum(w0 * theta_s) / sum(w0)
    Q_stat <- sum(w0 * ((theta_s - mu0)^2))
    denom_dl <- sum(w0) - sum(w0^2) / sum(w0)
    tau2_dl <- if (denom_dl > 1e-12) max(0.0, (Q_stat - (S_k - 1.0)) / denom_dl) else 0.0

    w_anc <- if ("n_trials" %in% names(st) && all(is.finite(st$n_trials)) && all(st$n_trials > 0)) {
      as.numeric(st$n_trials)
    } else {
      rep(1.0, S_k)
    }
    w_anc <- w_anc / sum(w_anc)
    mu_re <- sum(w_anc * theta_s)
    var_hksj <- sum((w_anc^2) * ((theta_s - mu_re)^2)) * (S_k / (S_k - 1.0))
    se_re <- sqrt(max(var_hksj, .Machine$double.eps))
    df_re <- S_k - 1.0
    t_re <- mu_re / se_re
    p_re <- 2.0 * stats::pt(abs(t_re), df = df_re, lower.tail = FALSE)
    t_crit <- stats::qt(0.975, df = df_re)
    grp_rows[[ti]] <- data.table::data.table(
      term = tm,
      estimate = mu_re,
      std.error = se_re,
      statistic = t_re,
      df = df_re,
      p.value = p_re,
      conf.low = mu_re - t_crit * se_re,
      conf.high = mu_re + t_crit * se_re,
      tau2_between = tau2_dl,
      mean_unweighted = mean(theta_raw),
      median_subj = stats::median(theta_raw),
      n_subjects = S_k
    )
  }

  grp_dt <- data.table::rbindlist(grp_rows, fill = TRUE)
  if (nrow(grp_dt) > 0L) {
    grp_dt[, `:=`(
      seed_region = seed_region,
      target_region = target_region,
      seed_lag = seed_lag,
      target_lag = target_lag,
      method = method
    )]
  }

  structure(
    list(
      group = grp_dt,
      subjects = subj_dt,
      panel_repair = if (isTRUE(return_augmented_panel)) data.table::rbindlist(aug_panels, fill = TRUE) else NULL,
      settings = list(
        seed_region = seed_region,
        target_region = target_region,
        seed_lag = seed_lag,
        target_lag = target_lag,
        covariates = covariates,
        test_ppi = test_ppi,
        control_target_seed_lag = control_target_seed_lag,
        orthogonalize_seed_x = orthogonalize_seed_x,
        method = method,
        use_whitened = use_whitened
      )
    ),
    class = "medusa_connectivity_repair"
  )
}

#' Between-subject REPAIR structural regression on MEDuSA mixed models
#'
#' Connects fitted \code{lme4::lmer()} (or \code{\link{mixed_by}} model objects)
#' from MEDuSA Stage 4 directly to \code{repair::repair()} (Vig & Hallquist,
#' 2026) to estimate measurement-error-corrected between-subject regressions
#' relating subject-specific neural random intercepts/slopes across brain
#' regions or to external behavioral/clinical variables.
#'
#' @param formula Structural regression formula using \code{repair::re()} markers
#'   (e.g., \code{re("x_whitened", model = "target") ~ re("x_whitened", model = "seed") + age}).
#' @param models Named list of fitted \code{lme4::lmer} or \code{nlme::lme}
#'   models (or a single fitted model).
#' @param data Subject-level \code{data.frame} with one row per subject ID.
#' @param id Character scalar naming the subject ID column in \code{data}
#'   (default \code{"id"}).
#' @param regularization Regularization policy passed to \code{repair::repair()}:
#'   \code{"step_up"} (default; adaptive M-MOMS), \code{"minimum"},
#'   \code{"none"}, or \code{"custom"}.
#' @param cross_model_error Cross-model prediction-error covariance policy when
#'   \code{models} has multiple entries (defaults to \code{"independent"} when
#'   \code{length(models) > 1L}).
#' @param medusa_mlm Optional alias for \code{models} (e.g., a single fitted
#'   \code{lme4::lmer} model).
#' @param subject_data Optional alias for \code{data}.
#' @param id_col Optional alias for \code{id}.
#' @param ... Additional arguments passed to \code{repair::repair()}.
#' @return A list with \code{coefficients} (tidy \code{data.frame} from
#'   \code{repair::repair_coefficients()}), \code{reliability} (from
#'   \code{repair::repair_reliability()}), \code{predictions} (likelihood-only
#'   random-effect predictions and error variances), and \code{fit} (the S7
#'   \code{class_repair_mod} object).
#' @export
medusa_between_subject_repair <- function(
  formula,
  models = NULL,
  data = NULL,
  id = "id",
  regularization = c("step_up", "minimum", "none", "custom"),
  cross_model_error = NULL,
  medusa_mlm = NULL,
  subject_data = NULL,
  id_col = NULL,
  ...
) {
  if (!requireNamespace("repair", quietly = TRUE)) {
    stop("The `repair` package is required for `medusa_between_subject_repair()`.")
  }
  if (is.null(models) && !is.null(medusa_mlm)) models <- medusa_mlm
  if (is.null(data) && !is.null(subject_data)) data <- subject_data
  if (!is.null(id_col)) id <- id_col
  if (is.null(models)) stop("`models` (or `medusa_mlm`) must be provided.")
  if (is.null(data)) stop("`data` (or `subject_data`) must be provided.")
  regularization <- match.arg(regularization)
  if (is.null(cross_model_error) && is.list(models) && length(models) > 1L) {
    cross_model_error <- "independent"
  }
  r_args <- list(
    formula = formula,
    models = models,
    data = as.data.frame(data),
    id = id,
    regularization = regularization
  )
  if (!is.null(cross_model_error)) r_args$cross_model_error <- cross_model_error
  extra_dots <- list(...)
  if (length(extra_dots) > 0L) r_args <- c(r_args, extra_dots)
  r_mod <- do.call(repair::repair, r_args)
  list(
    coefficients = repair::repair_coefficients(r_mod),
    reliability = repair::repair_reliability(r_mod),
    predictions = repair::repair_predictions(r_mod),
    fit = r_mod
  )
}

# Build canonical reference cross-lag resolution matrix C_sym and operator error
# correlation R_op when only a lag grid (and not a decon_epoch_regression object)
# is supplied to Stage-4 multi-lag inference.
medusa_reference_lag_operators <- function(
  lags,
  tr = 1.0,
  readout_kernel = "triangle",
  readout_width = 1.5,
  lambda = 0.8,
  rho = 0.30,
  penalty = "sobolev",
  gain_floor = 0.15
) {
  lags <- as.numeric(lags)
  n_lags <- length(lags)
  lag_names <- as.character(lags)
  if (n_lags < 2L) {
    I1 <- diag(1, max(1L, n_lags))
    if (n_lags == 1L) dimnames(I1) <- list(lag_names, lag_names)
    return(list(C_sym = I1, R_op = I1, d_lag = rep(1, n_lags)))
  }
  tr_eff <- if (is.finite(tr) && tr > 0.05) tr else 1.0
  w_eff <- if (!is.null(readout_width) && is.finite(readout_width[1L]) && readout_width[1L] > 0) {
    readout_width[1L]
  } else {
    1.5 * tr_eff
  }
  hrf_k <- hrf_library(tr = tr_eff, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6.0)[, 1L, drop = FALSE]
  span_lags <- max(lags) - min(lags)
  n_ref <- max(96L, as.integer(ceiling((span_lags + 80.0) / tr_eff)))
  rows_ref <- seq_len(n_ref)
  z_ref <- cbind(1, seq(-1, 1, length.out = n_ref))
  settings_ref <- list(
    lambda = lambda,
    scaling = "noise_sd",
    penalty = penalty,
    backend = "cpp",
    tail_padding = 0L,
    ridge_floor = 1e-8,
    nuisance_mode = "joint",
    observation_mask = "all"
  )
  op_ref <- hrf_inverse_operator(n_ref, hrf_k[, 1L], rows_ref, rho, z_ref, settings_ref)
  t_start <- max(12.0, abs(min(lags)) + 6.0)
  t_end <- max(t_start + 32.0, (n_ref - 16L) * tr_eff - max(0, max(lags)))
  onsets_ref <- seq(t_start, t_end, by = max(6.3, 4.3 * tr_eff))
  if (length(onsets_ref) < 4L) {
    onsets_ref <- seq(t_start, t_start + 30.0 * tr_eff, length.out = 6L)
  }
  m_ref <- length(onsets_ref)
  max_run_time <- (n_ref - 1L) * tr_eff

  S_lags_ref <- lapply(lags, function(lv) {
    build_epoch_readout_matrix(
      onsets_ref, lv, n_ref, tr_eff, 0,
      readout_kernel = readout_kernel, readout_width = w_eff
    )
  })
  ok_list <- lapply(S_lags_ref, function(S_l) rowSums(!is.finite(S_l)) == 0L)
  ok_common <- Reduce("&", ok_list)
  if (sum(ok_common) < 2L) {
    dist_mat <- abs(outer(lags, lags, "-"))
    C_fb <- exp(-0.5 * (dist_mat / (1.35 * tr_eff))^2)
    R_fb <- exp(-0.5 * (dist_mat / (2.10 * tr_eff))^2)
    dimnames(C_fb) <- dimnames(R_fb) <- list(lag_names, lag_names)
    return(list(C_sym = C_fb, R_op = R_fb, d_lag = rep(0.28, n_lags)))
  }
  m_c <- sum(ok_common)
  DG_lags_ref <- lapply(lags, function(lv) {
    sh_onsets <- pmin(max_run_time, pmax(0, onsets_ref + lv))
    G_u <- hrf_event_design_cpp(n_ref, tr_eff, 0, sh_onsets, rep(0, m_ref), seq_len(m_ref), m_ref, hrf_k)[[1L]]
    G_c <- sweep(G_u, 2L, colMeans(G_u), "-")
    op_ref %*% G_c
  })
  K_ref <- matrix(0, n_lags, n_lags)
  U_lag_mat <- matrix(0, n_ref * m_c, n_lags)
  for (l1 in seq_len(n_lags)) {
    S_l1 <- S_lags_ref[[l1]][ok_common, , drop = FALSE]
    for (l2 in seq_len(n_lags)) {
      DG_l2 <- DG_lags_ref[[l2]][, ok_common, drop = FALSE]
      k_diag <- rowSums(S_l1 * t(DG_l2))
      K_ref[l1, l2] <- mean(k_diag[is.finite(k_diag)])
    }
    B_l1 <- S_l1 %*% op_ref
    U_l1 <- compute_ar1_operator_factor(B_l1, rows_ref, rho, center_obs = TRUE)
    U_cent <- sweep(U_l1, 1L, rowMeans(U_l1), "-")
    U_lag_mat[, l1] <- as.numeric(U_cent)
  }
  d_ref <- pmax(diag(K_ref), gain_floor)
  D_inv <- diag(1 / sqrt(d_ref), n_lags)
  C_lag <- D_inv %*% K_ref %*% D_inv
  C_sym <- 0.5 * (C_lag + t(C_lag))
  V_lag <- crossprod(U_lag_mat) / max(1L, m_c - 1L)
  d_v <- sqrt(pmax(diag(V_lag), .Machine$double.eps))
  R_op <- V_lag / outer(d_v, d_v)
  R_op[!is.finite(R_op)] <- 0
  diag(R_op) <- 1
  dimnames(C_sym) <- dimnames(R_op) <- list(lag_names, lag_names)
  list(C_sym = C_sym, R_op = R_op, d_lag = d_ref)
}

# Build lobe-normalized spectral Lavrentiev-Tikhonov cross-lag de-blurring operator W_deb
medusa_build_deblur_operator <- function(
  C_sym,
  lags,
  tr = 1.0,
  alpha1 = 0.025,
  alpha0_sq = 0.0012,
  sigma_lobe = 1.05
) {
  n_lags <- length(lags)
  if (n_lags < 3L || is.null(C_sym) || !all(dim(C_sym) == c(n_lags, n_lags))) {
    return(diag(1, n_lags))
  }
  tr_eff <- if (is.finite(tr) && tr > 0.05) tr else 1.0
  dt_lag <- if (n_lags > 1L) stats::median(diff(sort(unique(lags)))) else tr_eff
  density_scale <- max(1.0, tr_eff / max(dt_lag, 0.1 * tr_eff))
  C_sym <- 0.5 * (C_sym + t(C_sym))
  eig_C <- eigen(C_sym, symmetric = TRUE)
  mu_sc <- pmax(eig_C$values, 0) / density_scale
  f_spec <- (mu_sc / (mu_sc^2 + alpha1 * mu_sc + alpha0_sq)) / density_scale
  W_raw <- eig_C$vectors %*% (f_spec * t(eig_C$vectors))
  R_raw <- W_raw %*% C_sym
  lag_dist <- abs(outer(lags, lags, "-"))
  L_lobe <- exp(-0.5 * (lag_dist / max(sigma_lobe, 0.25))^2)
  g_lobe <- pmax(rowSums(R_raw * L_lobe), 0.20)
  sweep(W_raw, 1L, g_lobe, "/")
}

# Joint cross-lag spectral de-blurring + correlated-SURE 2nd-order boosted P-spline smoother
medusa_smooth_deblur_trajectory <- function(
  beta_vec,
  Sigma_beta,
  lags,
  C_sym = NULL,
  tr = 1.0,
  demix_within_trial = FALSE,
  within_trial_ridge = 0.15,
  recalibrate_wt_lobe = FALSE,
  trajectory_mode = c("auto", "deblur_pspline", "pspline", "deblur", "none"),
  alpha1 = 0.025,
  alpha0_sq = 0.0012,
  sigma_lobe = 1.05
) {
  trajectory_mode <- match.arg(trajectory_mode)
  L <- length(beta_vec)
  if (L < 3L || identical(trajectory_mode, "none")) {
    se_out <- sqrt(pmax(diag(Sigma_beta), .Machine$double.eps))
    return(list(
      beta = as.numeric(beta_vec),
      se = se_out,
      Sigma = Sigma_beta,
      M_op = diag(1, L),
      eta_opt = 0,
      df_eff = L
    ))
  }
  do_deblur <- if (identical(trajectory_mode, "auto")) {
    !isTRUE(demix_within_trial) && !is.null(C_sym)
  } else {
    trajectory_mode %in% c("deblur_pspline", "deblur")
  }
  do_pspline <- trajectory_mode %in% c("auto", "deblur_pspline", "pspline")

  if (do_deblur && !is.null(C_sym)) {
    W_pre <- medusa_build_deblur_operator(
      C_sym = C_sym,
      lags = lags,
      tr = tr,
      alpha1 = alpha1,
      alpha0_sq = alpha0_sq,
      sigma_lobe = sigma_lobe
    )
  } else if (isTRUE(demix_within_trial) && isTRUE(recalibrate_wt_lobe) && !is.null(C_sym)) {
    W_1 <- solve(crossprod(C_sym) + within_trial_ridge * diag(L), t(C_sym))
    W_sym <- W_1 + 0.6 * W_1 %*% (diag(L) - C_sym %*% W_1)
    R_mat <- W_sym %*% C_sym
    diag_base <- pmax(diag(R_mat), 0.15)
    L_lobe <- exp(-0.5 * (outer(lags, lags, "-") / max(sigma_lobe, 0.25))^2)
    g_lobe_wt <- pmax(rowSums(R_mat * L_lobe), 0.20)
    W_pre <- diag(diag_base / g_lobe_wt, L)
  } else {
    W_pre <- diag(1, L)
  }

  b_pre <- as.numeric(W_pre %*% beta_vec)
  Sigma_pre <- W_pre %*% Sigma_beta %*% t(W_pre)
  Sigma_pre <- 0.5 * (Sigma_pre + t(Sigma_pre))

  if (!do_pspline || L < 4L) {
    se_out <- sqrt(pmax(diag(Sigma_pre), .Machine$double.eps))
    return(list(
      beta = b_pre,
      se = se_out,
      Sigma = Sigma_pre,
      M_op = W_pre,
      eta_opt = 0,
      df_eff = sum(diag(W_pre))
    ))
  }

  # Second-order difference penalty on potentially unequally spaced lags
  D2 <- matrix(0, L - 2L, L)
  for (i in seq_len(L - 2L)) {
    dt1 <- max(abs(lags[i + 1L] - lags[i]), 1e-4)
    dt2 <- max(abs(lags[i + 2L] - lags[i + 1L]), 1e-4)
    dt_bar <- 0.5 * (dt1 + dt2)
    D2[i, i:(i + 2L)] <- c(1 / (dt1 * dt_bar), -(1 / dt1 + 1 / dt2) / dt_bar, 1 / (dt2 * dt_bar))
  }
  P2 <- crossprod(D2)

  s2_diag <- pmax(diag(Sigma_pre), 1e-8)
  w_prec <- 1 / s2_diag
  w_prec <- w_prec / mean(w_prec)
  W_p <- diag(w_prec, L)

  eta_max <- if (do_deblur) 0.08 else 1.25
  eta_min <- if (do_deblur) 1e-4 else 2e-3
  eta_grid <- c(0, exp(seq(log(eta_min), log(eta_max), length.out = 28L)))
  IL <- diag(1, L)

  best_sure <- Inf
  best_S <- IL
  best_eta <- 0
  for (eta in eta_grid) {
    if (eta <= 0) {
      S_eta <- IL
    } else {
      S0 <- tryCatch(solve(W_p + eta * P2, W_p), error = function(e) IL)
      # Twicing / 2nd-order boosted P-spline reduces quadratic peak curvature bias from O(eta) to O(eta^2)
      S_eta <- S0 %*% (2 * IL - S0)
    }
    ImS <- IL - S_eta
    diff_b <- as.numeric(ImS %*% b_pre)
    rss_w <- sum(w_prec * (diff_b^2))
    # Exact generalized Stein's Unbiased Risk Estimate (SURE) under correlated cross-lag covariance Sigma_pre
    tr_1 <- sum(rowSums((W_p %*% ImS) * (ImS %*% Sigma_pre)))
    tr_2 <- sum(rowSums((W_p %*% S_eta) * (S_eta %*% Sigma_pre)))
    sure_val <- rss_w - tr_1 + tr_2
    if (is.finite(sure_val) && sure_val < best_sure) {
      best_sure <- sure_val
      best_S <- S_eta
      best_eta <- eta
    }
  }

  M_total <- best_S %*% W_pre
  b_post <- as.numeric(best_S %*% b_pre)
  Sigma_post <- best_S %*% Sigma_pre %*% t(best_S)
  Sigma_post <- 0.5 * (Sigma_post + t(Sigma_post))
  se_post <- sqrt(pmax(diag(Sigma_post), .Machine$double.eps))
  list(
    beta = b_post,
    se = se_post,
    Sigma = Sigma_post,
    M_op = M_total,
    eta_opt = best_eta,
    df_eff = sum(diag(best_S))
  )
}

# Vectorized 1D Threshold-Free Cluster Enhancement (Smith & Nichols, 2009) across a matrix
# of lag t-statistics T_mat (L x B). Uses a fixed integration step dh so observed and
# permutation/Monte-Carlo null columns are evaluated on the exact same threshold grid.
medusa_tfce_1d_matrix <- function(T_mat, dt = 1.0, E = 0.5, H = 2.0, dh = 0.10, two_sided = TRUE) {
  if (is.null(dim(T_mat))) T_mat <- matrix(T_mat, ncol = 1L)
  L <- nrow(T_mat)
  B <- ncol(T_mat)
  tfce_side <- function(M_pos) {
    max_val <- max(M_pos, na.rm = TRUE)
    if (!is.finite(max_val) || max_val < dh) return(matrix(0, L, B))
    h_max <- min(max_val, 14.0)
    h_grid <- seq(dh, h_max, by = dh)
    acc <- matrix(0, L, B)
    for (h in h_grid) {
      A <- M_pos >= h
      if (!any(A)) break
      # Forward accumulation of contiguous supra-threshold run lengths along lags
      F_len <- matrix(0L, L, B)
      F_len[1L, ] <- as.integer(A[1L, ])
      if (L >= 2L) {
        for (l in 2L:L) {
          al <- as.integer(A[l, ])
          F_len[l, ] <- al * (F_len[l - 1L, ] + 1L)
        }
        # Backward propagation so every lag in a contiguous cluster receives the full cluster extent
        for (l in (L - 1L):1L) {
          al <- as.integer(A[l, ])
          F_len[l, ] <- al * pmax(F_len[l, ], F_len[l + 1L, ])
        }
      }
      acc <- acc + ((F_len * dt)^E) * (h^H) * dh
    }
    acc
  }
  if (isTRUE(two_sided)) {
    tfce_pos <- tfce_side(pmax(T_mat, 0))
    tfce_neg <- tfce_side(pmax(-T_mat, 0))
    pmax(tfce_pos, tfce_neg)
  } else {
    tfce_side(pmax(T_mat, 0))
  }
}

# Continuous sub-grid natural cubic spline latency, center-of-mass, and FWHM with
# parametric bootstrap SEs and confidence intervals from the joint L x L covariance
medusa_spline_latency <- function(
  beta_vec,
  Sigma_beta,
  lags,
  fine_dt = 0.01,
  n_boot = 400L,
  conf_level = 0.95,
  polarity = c("positive", "negative", "abs")
) {
  polarity <- match.arg(polarity)
  L <- length(lags)
  if (L < 2L || any(!is.finite(beta_vec))) {
    pk <- if (L >= 1L && is.finite(beta_vec[1L])) lags[1L] else NA_real_
    return(list(
      peak_lag = pk,
      peak_lag_grid = pk,
      peak_lag_se = NA_real_,
      peak_lag_conf.low = NA_real_,
      peak_lag_conf.high = NA_real_,
      peak_val = if (L >= 1L) beta_vec[1L] else NA_real_,
      peak_val_se = NA_real_,
      com_lag = pk,
      com_lag_se = NA_real_,
      fwhm = NA_real_,
      fwhm_se = NA_real_
    ))
  }
  tau_fine <- seq(min(lags), max(lags), by = fine_dt)
  G <- length(tau_fine)
  # Construct linear natural cubic spline interpolation basis M_fine (G x L) once
  IL <- diag(1, L)
  M_fine <- vapply(seq_len(L), function(j) {
    stats::spline(x = lags, y = IL[, j], xout = tau_fine, method = "natural")$y
  }, numeric(G))
  if (is.null(dim(M_fine))) M_fine <- matrix(M_fine, nrow = G, ncol = L)

  extract_metrics <- function(f_curve) {
    s_curve <- if (identical(polarity, "positive")) {
      f_curve
    } else if (identical(polarity, "negative")) {
      -f_curve
    } else {
      abs(f_curve)
    }
    idx_pk <- which.max(s_curve)
    t_pk <- tau_fine[idx_pk]
    v_pk <- f_curve[idx_pk]
    s_pk <- s_curve[idx_pk]
    pos_w <- pmax(s_curve, 0)^2
    t_com <- if (sum(pos_w) > 1e-12) sum(tau_fine * pos_w) / sum(pos_w) else t_pk
    if (is.finite(s_pk) && s_pk > 1e-8) {
      above_half <- s_curve >= (0.5 * s_pk)
      # Contiguous lobe containing idx_pk
      left_idx <- idx_pk
      while (left_idx > 1L && above_half[left_idx - 1L]) left_idx <- left_idx - 1L
      right_idx <- idx_pk
      while (right_idx < G && above_half[right_idx + 1L]) right_idx <- right_idx + 1L
      fwhm_val <- tau_fine[right_idx] - tau_fine[left_idx]
    } else {
      fwhm_val <- NA_real_
    }
    c(peak_lag = t_pk, peak_val = v_pk, com_lag = t_com, fwhm = fwhm_val)
  }

  f_hat <- as.numeric(M_fine %*% beta_vec)
  obs_m <- extract_metrics(f_hat)
  grid_s <- if (identical(polarity, "positive")) {
    beta_vec
  } else if (identical(polarity, "negative")) {
    -beta_vec
  } else {
    abs(beta_vec)
  }
  peak_lag_grid <- lags[which.max(grid_s)]

  if (n_boot <= 0L || any(!is.finite(Sigma_beta))) {
    return(list(
      peak_lag = unname(obs_m["peak_lag"]),
      peak_lag_grid = peak_lag_grid,
      peak_lag_se = NA_real_,
      peak_lag_conf.low = NA_real_,
      peak_lag_conf.high = NA_real_,
      peak_val = unname(obs_m["peak_val"]),
      peak_val_se = NA_real_,
      com_lag = unname(obs_m["com_lag"]),
      com_lag_se = NA_real_,
      fwhm = unname(obs_m["fwhm"]),
      fwhm_se = NA_real_
    ))
  }

  Sigma_sym <- 0.5 * (Sigma_beta + t(Sigma_beta))
  eig_S <- eigen(Sigma_sym, symmetric = TRUE)
  L_half <- eig_S$vectors %*% (sqrt(pmax(eig_S$values, 0)) * t(eig_S$vectors))
  Z_draw <- matrix(stats::rnorm(L * n_boot), nrow = L, ncol = n_boot)
  B_draw <- beta_vec + L_half %*% Z_draw
  F_draw <- M_fine %*% B_draw
  boot_mat <- apply(F_draw, 2L, extract_metrics)

  alpha_tail <- (1.0 - conf_level) / 2.0
  pk_draws <- boot_mat["peak_lag", ]
  val_draws <- boot_mat["peak_val", ]
  com_draws <- boot_mat["com_lag", ]
  fwhm_draws <- boot_mat["fwhm", ]
  ci_pk <- stats::quantile(pk_draws, probs = c(alpha_tail, 1.0 - alpha_tail), na.rm = TRUE, names = FALSE)

  list(
    peak_lag = unname(obs_m["peak_lag"]),
    peak_lag_grid = peak_lag_grid,
    peak_lag_se = stats::sd(pk_draws, na.rm = TRUE),
    peak_lag_conf.low = ci_pk[1L],
    peak_lag_conf.high = ci_pk[2L],
    peak_val = unname(obs_m["peak_val"]),
    peak_val_se = stats::sd(val_draws, na.rm = TRUE),
    com_lag = unname(obs_m["com_lag"]),
    com_lag_se = stats::sd(com_draws, na.rm = TRUE),
    fwhm = unname(obs_m["fwhm"]),
    fwhm_se = stats::sd(fwhm_draws, na.rm = TRUE)
  )
}

#' Joint Multi-Lag Trajectory Modeling and Operator-Informed Inference for MEDuSA
#'
#' Replaces independent lag-by-lag post-hoc testing in MEDuSA Stage 4 with a
#' unified multi-lag estimation and inference framework that exploits the exact
#' Stage-2/3 cross-lag point-spread resolution matrix \eqn{C_{\mathrm{sym}}} and
#' cross-lag operator error correlation \eqn{R_{\mathrm{op}}}:
#' \enumerate{
#'   \item \strong{Lobe-Normalized Spectral Lavrentiev--Tikhonov Cross-Lag De-Blurring
#'     + Correlated-SURE 2nd-Order Boosted P-Spline Smoothing}: sharpens multiphasic
#'     and ramped within-trial neural waveforms without boundary distortion or
#'     interior gain inflation.
#'   \item \strong{Multivariate Empirical Bayes (EB) Subject Trajectory Shrinkage}:
#'     borrows strength across all \eqn{L} lags simultaneously via the between-subject
#'     trajectory covariance \eqn{\hat{\Sigma}_b} and subject-specific cross-lag
#'     error covariance \eqn{\hat{\Sigma}_{e,s}}.
#'   \item \strong{Operator-Informed Simultaneous Cross-Lag Inference}: calibrates
#'     epoch-wide family-wise error rate (FWER) via exact correlated multivariate-\eqn{t}
#'     Max-\eqn{|T|} simultaneous confidence bands (\code{p_maxT}, \code{conf.low.sim},
#'     \code{conf.high.sim}) and 1D Threshold-Free Cluster Enhancement (\code{tfce_stat},
#'     \code{p_tfce}), avoiding the severe conservatism of Bonferroni/BY-FDR under
#'     strong positive cross-lag correlation.
#'   \item \strong{Continuous Sub-Grid Spline Latency & Shape Estimation}: estimates
#'     continuous peak latency (\code{peak_lag}), center-of-mass latency (\code{com_lag}),
#'     and full-width at half-maximum (\code{fwhm}) with parametric-bootstrap SEs and
#'     confidence intervals.
#' }
#'
#' @param data Either: (a) a \code{\link{decon_epoch_regression}} result list (or a
#'   list of per-subject \code{decon_epoch_regression} results), (b) a trial-level
#'   MEDuSA panel \code{data.frame}/\code{data.table} (single-subject or multi-subject),
#'   or (c) a lag-level coefficient table with columns \code{lag} (or \code{evt_time}),
#'   \code{estimate} (or \code{beta_corrected}), and \code{std.error} (or \code{se_corrected}).
#' @param covariates Optional character vector of covariate names when \code{data} is a
#'   trial-level panel. If \code{NULL}, inferred automatically from panel columns.
#' @param lag_resolution Optional named list (by region/parcel/voxel) of \eqn{L \times L}
#'   cross-lag resolution matrices \eqn{C_{\mathrm{sym}}}, or a single \eqn{L \times L} matrix.
#' @param lag_error_cor Optional named list (by region/parcel/voxel) of \eqn{L \times L}
#'   operator cross-lag error correlation matrices \eqn{R_{\mathrm{op}}}, or a single matrix.
#' @param trajectory_mode Trajectory estimation mode: \code{"auto"} (default; applies
#'   spectral cross-lag de-blurring when \code{!demix_within_trial} followed by
#'   correlated-SURE boosted P-spline smoothing), \code{"deblur_pspline"},
#'   \code{"pspline"}, \code{"deblur"}, or \code{"none"}.
#' @param use_whitened Logical; if \code{TRUE} (default) and whitened columns
#'   (\code{decon_demixed_whitened}, \code{*_whitened}) are present in a trial panel,
#'   uses them for within-subject lag regressions.
#' @param tr Acquisition TR in seconds (default \code{1.0}), used when constructing
#'   reference cross-lag operators if \code{lag_resolution} is not attached.
#' @param readout_kernel Readout kernel (\code{"triangle"}, \code{"boxcar"},
#'   \code{"gaussian"}, \code{"linear"}) for reference operator construction.
#' @param readout_width Readout kernel width in seconds (default \code{1.5}).
#' @param decon_lambda Reference deconvolution penalty \eqn{\lambda} (default \code{0.8}).
#' @param rho_ar1 Reference AR(1) parameter (default \code{0.30}).
#' @param demix_within_trial Logical indicating whether Stage-3 within-trial cross-lag
#'   de-mixing was already applied to the input panel/coefficients (default \code{FALSE}).
#' @param recalibrate_wt_lobe Logical; if \code{TRUE} and \code{demix_within_trial = TRUE},
#'   recalibrates interior vs. boundary lobe gain on \code{W_lag}.
#' @param contrasts Optional named list of lag contrast specifications. Each entry can be
#'   either a numeric weight vector of length \eqn{L} or a list \code{list(early = c(t1, t2),
#'   late = c(t3, t4))} specifying two lag windows to contrast (\code{late - early}).
#' @param n_perm Number of correlated multivariate-\eqn{t} Monte Carlo draws for
#'   simultaneous Max-\eqn{|T|} and 1D TFCE calibration (default \code{4000L}).
#' @param n_boot_latency Number of parametric bootstrap draws for continuous spline
#'   latency SEs and confidence intervals (default \code{400L}).
#' @param fine_dt Sub-grid spline evaluation step in seconds (default \code{0.01}).
#' @param conf_level Simultaneous and pointwise confidence level (default \code{0.95}).
#' @param polarity Peak polarity for latency extraction (\code{"positive"},
#'   \code{"negative"}, or \code{"abs"}).
#' @param seed Optional random seed for reproducible Monte Carlo / bootstrap draws.
#' @return An object of class \code{"medusa_multilag"} containing:
#'   \describe{
#'     \item{\code{coefficients}}{\code{data.table} of group/unit lag estimates
#'       (\code{estimate}, \code{std.error}, \code{estimate_multilag},
#'       \code{std.error_multilag}), pointwise \eqn{t}/\eqn{p}-values, FDR-adjusted
#'       \eqn{p}-values (\code{p_fdr_bh}, \code{p_fdr_by}), simultaneous Max-\eqn{|T|}
#'       \eqn{p}-values and confidence bands (\code{p_maxT}, \code{conf.low.sim},
#'       \code{conf.high.sim}), 1D TFCE statistics and family-wise \eqn{p}-values
#'       (\code{tfce_stat}, \code{p_tfce}), and epoch omnibus tests
#'       (\code{omnibus_stat}, \code{omnibus_p}).}
#'     \item{\code{latency}}{\code{data.table} of continuous sub-grid spline peak
#'       latency (\code{peak_lag}, \code{peak_lag_se}, \code{peak_lag_conf.low},
#'       \code{peak_lag_conf.high}), discrete grid peak (\code{peak_lag_grid}),
#'       center-of-mass latency (\code{com_lag}, \code{com_lag_se}), and FWHM
#'       (\code{fwhm}, \code{fwhm_se}) per unit and term.}
#'     \item{\code{contrasts}}{\code{data.table} of covariance-propagated multi-lag
#'       window contrasts (\code{estimate}, \code{std.error}, \code{statistic},
#'       \code{df}, \code{p.value}, \code{conf.low}, \code{conf.high}).}
#'     \item{\code{subjects}}{\code{data.table} of subject-specific lag trajectories
#'       (\code{estimate_raw}, \code{estimate_smooth}, and Multivariate Empirical
#'       Bayes \code{estimate_eb}) when multi-subject data are supplied.}
#'     \item{\code{subject_latency}}{\code{data.table} of subject-specific continuous
#'       spline peak and center-of-mass latencies extracted from \code{estimate_eb}.}
#'     \item{\code{covariance}}{Nested list of \eqn{L \times L} raw and multi-lag
#'       trajectory covariance and correlation matrices per unit and term.}
#'   }
#' @export
medusa_multilag_infer <- function(
  data,
  covariates = NULL,
  lag_resolution = NULL,
  lag_error_cor = NULL,
  trajectory_mode = c("auto", "deblur_pspline", "pspline", "deblur", "none"),
  use_whitened = TRUE,
  tr = 1.0,
  readout_kernel = "triangle",
  readout_width = 1.5,
  decon_lambda = 0.8,
  rho_ar1 = 0.30,
  demix_within_trial = FALSE,
  recalibrate_wt_lobe = FALSE,
  contrasts = NULL,
  n_perm = 4000L,
  n_boot_latency = 400L,
  fine_dt = 0.01,
  conf_level = 0.95,
  polarity = c("positive", "negative", "abs"),
  seed = NULL
) {
  trajectory_mode <- match.arg(trajectory_mode)
  polarity <- match.arg(polarity)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # 1. Unpack input object (decon_epoch_regression result, list of results, panel, or coef table)
  if (is.list(data) && !is.data.frame(data)) {
    if (!is.null(data$coefficients) && is.data.frame(data$coefficients)) {
      # Single decon_epoch_regression output
      if (is.null(lag_resolution)) lag_resolution <- data$lag_resolution
      if (is.null(lag_error_cor)) lag_error_cor <- data$lag_error_cor
      if (!is.null(data$settings$demix_within_trial)) {
        demix_within_trial <- isTRUE(data$settings$demix_within_trial)
      }
      if (!is.null(data$settings$readout_kernel)) readout_kernel <- data$settings$readout_kernel
      if (!is.null(data$settings$readout_width)) readout_width <- data$settings$readout_width
      if (!is.null(data$panel) && is.data.frame(data$panel) && nrow(data$panel) > 0L) {
        dt_in <- data.table::as.data.table(data$panel)
        input_kind <- "panel"
      } else {
        dt_in <- data.table::as.data.table(data$coefficients)
        input_kind <- "coef"
      }
    } else if (length(data) >= 1L && is.list(data[[1L]]) &&
               (!is.null(data[[1L]]$panel) || !is.null(data[[1L]]$coefficients))) {
      # List of per-subject decon_epoch_regression or fit_subject_modern_medusa outputs
      if (is.null(lag_resolution) && !is.null(data[[1L]]$lag_resolution)) {
        lag_resolution <- data[[1L]]$lag_resolution
      }
      if (is.null(lag_error_cor) && !is.null(data[[1L]]$lag_error_cor)) {
        lag_error_cor <- data[[1L]]$lag_error_cor
      }
      if (!is.null(data[[1L]]$settings$demix_within_trial)) {
        demix_within_trial <- isTRUE(data[[1L]]$settings$demix_within_trial)
      }
      has_panels <- all(vapply(data, function(x) !is.null(x$panel) && is.data.frame(x$panel), logical(1L)))
      if (has_panels) {
        sub_names <- names(data)
        if (is.null(sub_names)) sub_names <- sprintf("sub_%03d", seq_along(data))
        dt_in <- data.table::rbindlist(lapply(seq_along(data), function(i) {
          p_i <- data.table::as.data.table(data[[i]]$panel)
          if (!("id" %in% names(p_i)) && !("subject" %in% names(p_i))) {
            p_i[, id := sub_names[i]]
          }
          p_i
        }), fill = TRUE)
        input_kind <- "panel"
      } else {
        sub_names <- names(data)
        if (is.null(sub_names)) sub_names <- sprintf("sub_%03d", seq_along(data))
        dt_in <- data.table::rbindlist(lapply(seq_along(data), function(i) {
          c_i <- data.table::as.data.table(data[[i]]$coefficients)
          if (!("id" %in% names(c_i)) && !("subject" %in% names(c_i))) {
            c_i[, id := sub_names[i]]
          }
          c_i
        }), fill = TRUE)
        input_kind <- "coef"
      }
    } else {
      stop("Unsupported list input to `medusa_multilag_infer()`.")
    }
  } else if (is.data.frame(data)) {
    if (is.null(lag_resolution)) lag_resolution <- attr(data, "lag_resolution")
    if (is.null(lag_error_cor)) lag_error_cor <- attr(data, "lag_error_cor")
    dt_in <- data.table::as.data.table(data)
    if (any(c("decon_demixed", "decon_interp", "decon_raw", "decon_demixed_whitened") %in% names(dt_in))) {
      input_kind <- "panel"
    } else if (any(c("beta_corrected", "estimate", "beta_weighted", "beta_raw") %in% names(dt_in))) {
      input_kind <- "coef"
    } else {
      stop("Data frame supplied to `medusa_multilag_infer()` must be a MEDuSA trial panel or coefficient table.")
    }
  } else {
    stop("`data` must be a `decon_epoch_regression` result, list of results, or `data.frame`/`data.table`.")
  }

  # Identify unit column (region / parcel / voxel)
  unit_col <- intersect(c("region", "parcel", "voxel", "unit"), names(dt_in))[1L]
  if (is.na(unit_col)) {
    dt_in[, unit := "unit_1"]
    unit_col <- "unit"
  }
  lag_col <- intersect(c("evt_time", "lag"), names(dt_in))[1L]
  if (is.na(lag_col)) stop("Could not find `evt_time` or `lag` column in `data`.")
  subj_col <- intersect(c("id", "subject", "sub"), names(dt_in))[1L]

  lags <- sort(unique(as.numeric(dt_in[[lag_col]])))
  L <- length(lags)
  if (L < 2L) stop("`medusa_multilag_infer()` requires at least 2 lags.")
  dt_lag <- stats::median(diff(lags))

  # Construct or retrieve reference operators C_sym and R_op
  ref_ops <- NULL
  get_unit_operators <- function(u_id) {
    C_u <- if (is.matrix(lag_resolution) && all(dim(lag_resolution) == c(L, L))) {
      lag_resolution
    } else if (is.list(lag_resolution) && !is.null(lag_resolution[[as.character(u_id)]]) &&
               all(dim(lag_resolution[[as.character(u_id)]]) == c(L, L))) {
      lag_resolution[[as.character(u_id)]]
    } else {
      if (is.null(ref_ops)) {
        ref_ops <<- medusa_reference_lag_operators(
          lags = lags, tr = tr, readout_kernel = readout_kernel,
          readout_width = readout_width, lambda = decon_lambda, rho = rho_ar1
        )
      }
      ref_ops$C_sym
    }
    R_u <- if (is.matrix(lag_error_cor) && all(dim(lag_error_cor) == c(L, L))) {
      lag_error_cor
    } else if (is.list(lag_error_cor) && !is.null(lag_error_cor[[as.character(u_id)]]) &&
               all(dim(lag_error_cor[[as.character(u_id)]]) == c(L, L))) {
      lag_error_cor[[as.character(u_id)]]
    } else {
      if (is.null(ref_ops)) {
        ref_ops <<- medusa_reference_lag_operators(
          lags = lags, tr = tr, readout_kernel = readout_kernel,
          readout_width = readout_width, lambda = decon_lambda, rho = rho_ar1
        )
      }
      ref_ops$R_op
    }
    list(C_sym = C_u, R_op = R_u)
  }

  # 2. If input is a trial-level panel, extract subject-by-lag slope and SE matrices
  #    plus within-subject trial score cross-lag correlations
  if (identical(input_kind, "panel")) {
    has_w_y <- isTRUE(use_whitened) && ("decon_demixed_whitened" %in% names(dt_in)) &&
      any(is.finite(dt_in$decon_demixed_whitened))
    y_col <- if (has_w_y) {
      "decon_demixed_whitened"
    } else if ("decon_demixed" %in% names(dt_in) && any(is.finite(dt_in$decon_demixed))) {
      "decon_demixed"
    } else if ("decon_interp" %in% names(dt_in) && any(is.finite(dt_in$decon_interp))) {
      "decon_interp"
    } else {
      "decon_raw"
    }
    w_col <- if ("precision_weight" %in% names(dt_in)) "precision_weight" else NULL

    if (is.null(covariates)) {
      reserved_cols <- c(
        "id", "subject", "sub", "run", "trial", "event_index", "onset", "evt_time", "lag",
        "voxel", "parcel", "region", "unit", "decon_raw", "decon_demixed", "decon_eb",
        "decon_interp", "se_raw", "se_demixed", "reliability_trial", "precision_weight",
        "precision_weight_op", "precision_weight_raw", "error_ac1", "tau2_bio",
        "decon_demixed_whitened", "se_whitened", "decon_eb_whitened"
      )
      cand_covs <- setdiff(names(dt_in), reserved_cols)
      cand_covs <- cand_covs[!grepl("_whitened$", cand_covs)]
      cand_covs <- cand_covs[vapply(cand_covs, function(cn) is.numeric(dt_in[[cn]]), logical(1L))]
      covariates <- if (length(cand_covs) > 0L) cand_covs else "intercept"
    }
  } else {
    est_col <- intersect(c("beta_corrected", "estimate", "beta_weighted", "beta_raw"), names(dt_in))[1L]
    se_col <- intersect(c("se_corrected", "std.error", "se_weighted", "se_unweighted"), names(dt_in))[1L]
    if (is.na(est_col)) stop("Could not find coefficient estimate column in `data`.")
    if (!("term" %in% names(dt_in))) dt_in[, term := "x"]
    covariates <- unique(as.character(dt_in$term))
  }

  units <- unique(as.character(dt_in[[unit_col]]))
  has_multi_subj <- !is.na(subj_col) && length(unique(dt_in[[subj_col]])) >= 2L

  coef_out_list <- list()
  lat_out_list <- list()
  con_out_list <- list()
  subj_out_list <- list()
  subj_lat_list <- list()
  cov_out_list <- list()

  for (u_id in units) {
    dt_u <- dt_in[as.character(get(unit_col)) == u_id]
    ops_u <- get_unit_operators(u_id)
    C_u <- ops_u$C_sym
    R_op_u <- ops_u$R_op
    cov_out_list[[u_id]] <- list()

    for (tm in covariates) {
      # Extract B_mat (S x L), SE_mat (S x L), df_vec, and empirical residual/score correlation R_emp
      if (identical(input_kind, "panel")) {
        s_ids <- if (!is.na(subj_col)) unique(as.character(dt_u[[subj_col]])) else "sub_1"
        S_n <- length(s_ids)
        B_mat <- matrix(NA_real_, S_n, L, dimnames = list(s_ids, as.character(lags)))
        SE_mat <- matrix(NA_real_, S_n, L, dimnames = list(s_ids, as.character(lags)))
        df_sub <- numeric(S_n)
        score_cor_accum <- matrix(0, L, L)
        score_cor_n <- 0L

        tm_w <- paste0(tm, "_whitened")
        x_col_use <- if (has_w_y && (tm_w %in% names(dt_u))) tm_w else tm

        for (si in seq_len(S_n)) {
          dt_su <- if (!is.na(subj_col)) dt_u[as.character(get(subj_col)) == s_ids[si]] else dt_u
          # Matrix of trial-wise regression scores to capture within-subject cross-lag error correlation
          has_run_trial <- all(c("run", "trial") %in% names(dt_su))
          if (has_run_trial) {
            rt_keys <- unique(paste(dt_su$run, dt_su$trial, sep = "_"))
            score_mat <- matrix(NA_real_, length(rt_keys), L)
          } else {
            score_mat <- NULL
          }
          df_l_vec <- numeric(L)
          for (li in seq_len(L)) {
            dt_sul <- dt_su[abs(as.numeric(get(lag_col)) - lags[li]) < 1e-6]
            if (nrow(dt_sul) < 4L) next
            yv <- as.numeric(dt_sul[[y_col]])
            wv <- if (!is.null(w_col)) as.numeric(dt_sul[[w_col]]) else rep(1.0, length(yv))
            wv[!is.finite(wv) | wv <= 0] <- 1.0
            if (identical(tm, "intercept")) {
              ok_i <- is.finite(yv)
              if (sum(ok_i) < 3L) next
              w_ok <- wv[ok_i] / mean(wv[ok_i])
              y_ok <- yv[ok_i]
              b_hat <- sum(w_ok * y_ok) / sum(w_ok)
              res_i <- y_ok - b_hat
              df_i <- max(1L, sum(ok_i) - 1L)
              se_hat <- sqrt(sum((w_ok^2) * (res_i^2)) / (df_i * (sum(w_ok)^2)) * sum(ok_i))
              B_mat[si, li] <- b_hat
              SE_mat[si, li] <- se_hat
              df_l_vec[li] <- df_i
              if (!is.null(score_mat)) {
                k_ok <- paste(dt_sul$run[ok_i], dt_sul$trial[ok_i], sep = "_")
                score_mat[match(k_ok, rt_keys), li] <- w_ok * res_i
              }
            } else {
              xv <- as.numeric(dt_sul[[x_col_use]])
              ok_i <- is.finite(yv) & is.finite(xv)
              if (sum(ok_i) < 4L) next
              w_ok <- wv[ok_i] / mean(wv[ok_i])
              y_ok <- yv[ok_i]
              x_ok <- xv[ok_i]
              # Center within run if run column is present
              if ("run" %in% names(dt_sul)) {
                r_ok <- dt_sul$run[ok_i]
                for (rk in unique(r_ok)) {
                  idx_r <- which(r_ok == rk)
                  if (length(idx_r) >= 2L) {
                    wr <- w_ok[idx_r]
                    y_ok[idx_r] <- y_ok[idx_r] - sum(wr * y_ok[idx_r]) / sum(wr)
                    x_ok[idx_r] <- x_ok[idx_r] - sum(wr * x_ok[idx_r]) / sum(wr)
                  }
                }
              } else {
                y_ok <- y_ok - sum(w_ok * y_ok) / sum(w_ok)
                x_ok <- x_ok - sum(w_ok * x_ok) / sum(w_ok)
              }
              sxx <- sum(w_ok * (x_ok^2))
              if (sxx <= 1e-10) next
              b_hat <- sum(w_ok * x_ok * y_ok) / sxx
              res_i <- y_ok - x_ok * b_hat
              n_runs_i <- if ("run" %in% names(dt_sul)) length(unique(dt_sul$run[ok_i])) else 1L
              df_i <- max(1L, sum(ok_i) - 1L - n_runs_i)
              # White/sandwich variance combined with model SE for small-sample stability
              s2_ols <- sum(w_ok * (res_i^2)) / df_i
              v_model <- s2_ols / sxx
              v_hc <- sum((w_ok * x_ok * res_i)^2) / (sxx^2) * (sum(ok_i) / df_i)
              se_hat <- sqrt(max(0.5 * v_model + 0.5 * v_hc, .Machine$double.eps))
              B_mat[si, li] <- b_hat
              SE_mat[si, li] <- se_hat
              df_l_vec[li] <- df_i
              if (!is.null(score_mat)) {
                k_ok <- paste(dt_sul$run[ok_i], dt_sul$trial[ok_i], sep = "_")
                score_mat[match(k_ok, rt_keys), li] <- (w_ok * x_ok * res_i) / sxx
              }
            }
          }
          df_sub[si] <- stats::median(df_l_vec[df_l_vec > 0], na.rm = TRUE)
          if (!is.null(score_mat)) {
            ok_rows <- rowSums(!is.finite(score_mat)) == 0L
            if (sum(ok_rows) >= 5L) {
              R_sc <- suppressWarnings(stats::cor(score_mat[ok_rows, , drop = FALSE]))
              if (all(is.finite(R_sc))) {
                score_cor_accum <- score_cor_accum + R_sc
                score_cor_n <- score_cor_n + 1L
              }
            }
          }
        }
        R_within <- if (score_cor_n >= 1L) {
          0.65 * (score_cor_accum / score_cor_n) + 0.35 * R_op_u
        } else {
          R_op_u
        }
        diag(R_within) <- 1.0
      } else {
        # Coefficient table input
        dt_ut <- dt_u[as.character(term) == tm]
        if (nrow(dt_ut) == 0L) next
        s_ids <- if (!is.na(subj_col)) unique(as.character(dt_ut[[subj_col]])) else "sub_1"
        S_n <- length(s_ids)
        B_mat <- matrix(NA_real_, S_n, L, dimnames = list(s_ids, as.character(lags)))
        SE_mat <- matrix(NA_real_, S_n, L, dimnames = list(s_ids, as.character(lags)))
        df_sub <- rep(24, S_n)
        for (si in seq_len(S_n)) {
          dt_s <- if (!is.na(subj_col)) dt_ut[as.character(get(subj_col)) == s_ids[si]] else dt_ut
          ord <- match(lags, as.numeric(dt_s[[lag_col]]))
          B_mat[si, ] <- as.numeric(dt_s[[est_col]])[ord]
          if (!is.na(se_col)) {
            SE_mat[si, ] <- as.numeric(dt_s[[se_col]])[ord]
          } else {
            SE_mat[si, ] <- 0.10
          }
          if ("df" %in% names(dt_s)) {
            df_vals <- as.numeric(dt_s$df)[ord]
            if (any(is.finite(df_vals) & df_vals > 0)) {
              df_sub[si] <- stats::median(df_vals[is.finite(df_vals) & df_vals > 0])
            }
          }
        }
        R_within <- R_op_u
      }

      ok_subj <- rowSums(!is.finite(B_mat)) == 0L
      if (!any(ok_subj)) next
      B_mat <- B_mat[ok_subj, , drop = FALSE]
      SE_mat <- SE_mat[ok_subj, , drop = FALSE]
      df_sub <- df_sub[ok_subj]
      s_ids_ok <- rownames(B_mat)
      S_n <- nrow(B_mat)
      for (li in seq_len(L)) {
        bad_se <- !is.finite(SE_mat[, li]) | SE_mat[, li] <= 0
        if (any(bad_se)) {
          SE_mat[bad_se, li] <- if (S_n >= 2L) stats::sd(B_mat[, li], na.rm = TRUE) else 0.10
        }
      }

      # 3. Group / Unit Mean Trajectory & Full L x L Cross-Lag Covariance
      if (S_n >= 2L) {
        # Random-effects precision weights across subjects per lag
        b_raw_grp <- numeric(L)
        se_raw_grp <- numeric(L)
        W_subj <- matrix(1 / S_n, S_n, L)
        for (li in seq_len(L)) {
          yi <- B_mat[, li]
          s2i <- SE_mat[, li]^2
          v_tot <- stats::var(yi)
          tau2_b <- max(0, v_tot - mean(s2i))
          wi <- 1 / pmax(tau2_b + s2i, 1e-8)
          wi <- wi / sum(wi)
          W_subj[, li] <- wi
          b_raw_grp[li] <- sum(wi * yi)
          # Standard error of precision-weighted group mean
          res_y <- yi - b_raw_grp[li]
          v_emp <- sum((wi^2) * (res_y^2)) * (S_n / max(1L, S_n - 1L))
          v_mod <- sum((wi^2) * (tau2_b + s2i))
          se_raw_grp[li] <- sqrt(max(0.75 * v_emp + 0.25 * v_mod, .Machine$double.eps))
        }
        # Empirical between-subject cross-lag correlation shrunk toward operator correlation R_within
        B_cent <- sweep(B_mat, 2L, b_raw_grp, "-")
        R_emp_grp <- suppressWarnings(stats::cor(B_cent))
        if (any(!is.finite(R_emp_grp))) R_emp_grp <- R_within
        shrink_alpha <- min(0.75, max(0.20, (S_n - 2.0) / (S_n + 6.0)))
        R_grp <- shrink_alpha * R_emp_grp + (1.0 - shrink_alpha) * R_within
        R_grp <- 0.5 * (R_grp + t(R_grp))
        diag(R_grp) <- 1.0
        eig_R <- eigen(R_grp, symmetric = TRUE)
        R_grp <- eig_R$vectors %*% (pmax(eig_R$values, 0.02) * t(eig_R$vectors))
        d_rg <- sqrt(pmax(diag(R_grp), .Machine$double.eps))
        R_grp <- R_grp / outer(d_rg, d_rg)
        diag(R_grp) <- 1.0
        Sigma_raw_grp <- outer(se_raw_grp, se_raw_grp) * R_grp
        df_grp <- max(2, S_n - 1)
      } else {
        b_raw_grp <- as.numeric(B_mat[1L, ])
        se_raw_grp <- as.numeric(SE_mat[1L, ])
        R_grp <- R_within
        Sigma_raw_grp <- outer(se_raw_grp, se_raw_grp) * R_grp
        df_grp <- if (is.finite(df_sub[1L]) && df_sub[1L] >= 2) df_sub[1L] else 20
      }

      # 4. Joint Cross-Lag Spectral De-Blurring + Correlated-SURE Boosted P-Spline
      sm_grp <- medusa_smooth_deblur_trajectory(
        beta_vec = b_raw_grp,
        Sigma_beta = Sigma_raw_grp,
        lags = lags,
        C_sym = C_u,
        tr = tr,
        demix_within_trial = demix_within_trial,
        recalibrate_wt_lobe = recalibrate_wt_lobe,
        trajectory_mode = trajectory_mode
      )
      b_ml_grp <- sm_grp$beta
      se_ml_grp <- sm_grp$se
      Sigma_ml_grp <- sm_grp$Sigma

      # 5. Subject-Level Smooth and Multivariate Empirical Bayes (EB) Trajectories
      if (S_n >= 1L) {
        B_smooth_subj <- matrix(0, S_n, L)
        Sigma_e_subj <- vector("list", S_n)
        for (si in seq_len(S_n)) {
          Sigma_raw_si <- outer(SE_mat[si, ], SE_mat[si, ]) * R_within
          sm_si <- medusa_smooth_deblur_trajectory(
            beta_vec = B_mat[si, ],
            Sigma_beta = Sigma_raw_si,
            lags = lags,
            C_sym = C_u,
            tr = tr,
            demix_within_trial = demix_within_trial,
            recalibrate_wt_lobe = recalibrate_wt_lobe,
            trajectory_mode = trajectory_mode
          )
          B_smooth_subj[si, ] <- sm_si$beta
          Sigma_e_subj[[si]] <- sm_si$Sigma
        }

        if (S_n >= 2L) {
          Sigma_tot_sub <- stats::cov(B_smooth_subj)
          Sigma_e_bar <- Reduce("+", Sigma_e_subj) / S_n
          Diff_b <- 0.5 * ((Sigma_tot_sub - Sigma_e_bar) + t(Sigma_tot_sub - Sigma_e_bar))
          eig_b <- eigen(Diff_b, symmetric = TRUE)
          floor_b <- 0.05 * stats::median(diag(Sigma_e_bar))
          Sigma_b_hat <- eig_b$vectors %*% (pmax(eig_b$values, floor_b) * t(eig_b$vectors))

          B_eb_subj <- matrix(0, S_n, L)
          for (si in seq_len(S_n)) {
            S_tot_i <- Sigma_b_hat + Sigma_e_subj[[si]] + 1e-6 * diag(L)
            K_eb_i <- Sigma_b_hat %*% solve(S_tot_i)
            B_eb_subj[si, ] <- b_ml_grp + as.numeric(K_eb_i %*% (B_smooth_subj[si, ] - b_ml_grp))
          }
        } else {
          B_eb_subj <- B_smooth_subj
        }

        for (si in seq_len(S_n)) {
          df_s_rows <- data.table::data.table(
            id = s_ids_ok[si],
            unit = u_id,
            term = tm,
            lag = lags,
            estimate_raw = B_mat[si, ],
            std.error_raw = SE_mat[si, ],
            estimate_smooth = B_smooth_subj[si, ],
            estimate_eb = B_eb_subj[si, ]
          )
          df_s_rows[[unit_col]] <- u_id
          subj_out_list[[length(subj_out_list) + 1L]] <- df_s_rows

          lat_si <- medusa_spline_latency(
            beta_vec = B_eb_subj[si, ],
            Sigma_beta = Sigma_e_subj[[si]],
            lags = lags,
            fine_dt = fine_dt,
            n_boot = 0L,
            polarity = polarity
          )
          dt_lat_si <- data.table::data.table(
            id = s_ids_ok[si],
            unit = u_id,
            term = tm,
            peak_lag = lat_si$peak_lag,
            peak_lag_grid = lat_si$peak_lag_grid,
            peak_val = lat_si$peak_val,
            com_lag = lat_si$com_lag,
            fwhm = lat_si$fwhm
          )
          dt_lat_si[[unit_col]] <- u_id
          subj_lat_list[[length(subj_lat_list) + 1L]] <- dt_lat_si
        }
      }

      # 6. Operator-Informed Simultaneous Cross-Lag Inference (Max-|T|, 1D TFCE, & Omnibus)
      # Use raw lag-wise t-statistics coupled through the operator/empirical cross-lag
      # correlation R_grp for pointwise & cluster-aware hypothesis testing, while also
      # providing simultaneous bands for both raw and multi-lag trajectories.
      t_obs <- b_raw_grp / pmax(se_raw_grp, .Machine$double.eps)
      p_point <- 2.0 * stats::pt(abs(t_obs), df = df_grp, lower.tail = FALSE)
      p_bh <- stats::p.adjust(p_point, method = "BH")
      p_by <- stats::p.adjust(p_point, method = "BY")

      t_ml_obs <- b_ml_grp / pmax(se_ml_grp, .Machine$double.eps)
      p_ml_point <- 2.0 * stats::pt(abs(t_ml_obs), df = df_grp, lower.tail = FALSE)

      eig_Rg <- eigen(R_grp, symmetric = TRUE)
      L_Rg <- eig_Rg$vectors %*% (sqrt(pmax(eig_Rg$values, 0)) * t(eig_Rg$vectors))
      Z_null <- matrix(stats::rnorm(L * n_perm), nrow = L, ncol = n_perm)
      chi_draws <- sqrt(stats::rchisq(n_perm, df = df_grp) / df_grp)
      T_null <- sweep(L_Rg %*% Z_null, 2L, chi_draws, "/")

      max_t_null <- apply(abs(T_null), 2L, max)
      crit_maxT <- as.numeric(stats::quantile(max_t_null, probs = conf_level, na.rm = TRUE, names = FALSE))
      p_maxT <- vapply(abs(t_obs), function(tv) (1.0 + sum(max_t_null >= tv)) / (n_perm + 1.0), numeric(1L))

      # 1D Threshold-Free Cluster Enhancement (TFCE) calibrated under correlated T_null
      T_all <- cbind(t_obs, T_null)
      tfce_all <- medusa_tfce_1d_matrix(T_all, dt = dt_lag, E = 0.5, H = 2.0, dh = 0.10, two_sided = TRUE)
      tfce_obs <- tfce_all[, 1L]
      max_tfce_null <- apply(tfce_all[, -1L, drop = FALSE], 2L, max)
      p_tfce <- vapply(tfce_obs, function(sv) (1.0 + sum(max_tfce_null >= sv)) / (n_perm + 1.0), numeric(1L))

      # Regularized Mahalanobis / Hotelling epoch omnibus test against H0: beta(tau) == 0 for all tau
      eig_Sig <- eigen(Sigma_raw_grp, symmetric = TRUE)
      lam_floor <- 0.03 * max(eig_Sig$values, 1e-8)
      Sig_inv_reg <- eig_Sig$vectors %*% ((1.0 / pmax(eig_Sig$values, lam_floor)) * t(eig_Sig$vectors))
      q_obs <- sum(b_raw_grp * as.numeric(Sig_inv_reg %*% b_raw_grp)) / L
      B_null_draws <- sweep(T_null, 1L, se_raw_grp, "*")
      q_null <- colSums(B_null_draws * (Sig_inv_reg %*% B_null_draws)) / L
      omnibus_p <- (1.0 + sum(q_null >= q_obs)) / (n_perm + 1.0)

      t_crit_pt <- stats::qt(1.0 - (1.0 - conf_level) / 2.0, df = df_grp)
      dt_coef <- data.table::data.table(
        unit = u_id,
        term = tm,
        lag = lags,
        evt_time = lags,
        estimate = b_raw_grp,
        std.error = se_raw_grp,
        estimate_multilag = b_ml_grp,
        std.error_multilag = se_ml_grp,
        statistic = t_obs,
        statistic_multilag = t_ml_obs,
        df = df_grp,
        p.value = p_point,
        p.value_multilag = p_ml_point,
        p_fdr_bh = p_bh,
        p_fdr_by = p_by,
        p_maxT = p_maxT,
        tfce_stat = tfce_obs,
        p_tfce = p_tfce,
        crit_maxT = crit_maxT,
        conf.low = b_raw_grp - t_crit_pt * se_raw_grp,
        conf.high = b_raw_grp + t_crit_pt * se_raw_grp,
        conf.low.sim = b_raw_grp - crit_maxT * se_raw_grp,
        conf.high.sim = b_raw_grp + crit_maxT * se_raw_grp,
        conf.low.multilag.sim = b_ml_grp - crit_maxT * se_ml_grp,
        conf.high.multilag.sim = b_ml_grp + crit_maxT * se_ml_grp,
        omnibus_stat = q_obs,
        omnibus_p = omnibus_p,
        n_subjects = S_n
      )
      dt_coef[[unit_col]] <- u_id
      coef_out_list[[length(coef_out_list) + 1L]] <- dt_coef

      # 7. Continuous Sub-Grid Spline Peak Latency, Center-of-Mass, & FWHM
      lat_grp <- medusa_spline_latency(
        beta_vec = b_ml_grp,
        Sigma_beta = Sigma_ml_grp,
        lags = lags,
        fine_dt = fine_dt,
        n_boot = n_boot_latency,
        conf_level = conf_level,
        polarity = polarity
      )
      dt_lat <- data.table::data.table(
        unit = u_id,
        term = tm,
        peak_lag = lat_grp$peak_lag,
        peak_lag_grid = lat_grp$peak_lag_grid,
        peak_lag_se = lat_grp$peak_lag_se,
        peak_lag_conf.low = lat_grp$peak_lag_conf.low,
        peak_lag_conf.high = lat_grp$peak_lag_conf.high,
        peak_val = lat_grp$peak_val,
        peak_val_se = lat_grp$peak_val_se,
        com_lag = lat_grp$com_lag,
        com_lag_se = lat_grp$com_lag_se,
        fwhm = lat_grp$fwhm,
        fwhm_se = lat_grp$fwhm_se,
        omnibus_p = omnibus_p,
        min_p_maxT = min(p_maxT, na.rm = TRUE),
        min_p_tfce = min(p_tfce, na.rm = TRUE)
      )
      dt_lat[[unit_col]] <- u_id
      lat_out_list[[length(lat_out_list) + 1L]] <- dt_lat

      # 8. Covariance-Propagated Multi-Lag Contrasts
      con_specs <- contrasts
      if (is.null(con_specs) && L >= 4L) {
        mid_tau <- stats::median(lags)
        early_idx <- which(lags < mid_tau)
        late_idx <- which(lags > mid_tau)
        if (length(early_idx) >= 1L && length(late_idx) >= 1L) {
          c_default <- numeric(L)
          c_default[early_idx] <- -1.0 / length(early_idx)
          c_default[late_idx] <- 1.0 / length(late_idx)
          con_specs <- list(late_minus_early = c_default)
        }
      }
      if (!is.null(con_specs) && length(con_specs) > 0L) {
        c_names <- names(con_specs)
        if (is.null(c_names)) c_names <- paste0("contrast_", seq_along(con_specs))
        for (ci in seq_along(con_specs)) {
          cs <- con_specs[[ci]]
          if (is.numeric(cs) && length(cs) == L) {
            c_vec <- as.numeric(cs)
          } else if (is.list(cs) && !is.null(cs$early) && !is.null(cs$late)) {
            c_vec <- numeric(L)
            e_idx <- which(lags >= min(cs$early) - 1e-6 & lags <= max(cs$early) + 1e-6)
            l_idx <- which(lags >= min(cs$late) - 1e-6 & lags <= max(cs$late) + 1e-6)
            if (length(e_idx)) c_vec[e_idx] <- -1.0 / length(e_idx)
            if (length(l_idx)) c_vec[l_idx] <- 1.0 / length(l_idx)
          } else {
            next
          }
          est_c <- sum(c_vec * b_ml_grp)
          var_c <- as.numeric(crossprod(c_vec, Sigma_ml_grp %*% c_vec))
          se_c <- sqrt(max(var_c, .Machine$double.eps))
          t_c <- est_c / se_c
          p_c <- 2.0 * stats::pt(abs(t_c), df = df_grp, lower.tail = FALSE)
          dt_c <- data.table::data.table(
            unit = u_id,
            term = tm,
            contrast = c_names[ci],
            estimate = est_c,
            std.error = se_c,
            statistic = t_c,
            df = df_grp,
            p.value = p_c,
            conf.low = est_c - t_crit_pt * se_c,
            conf.high = est_c + t_crit_pt * se_c
          )
          dt_c[[unit_col]] <- u_id
          con_out_list[[length(con_out_list) + 1L]] <- dt_c
        }
      }

      cov_out_list[[u_id]][[tm]] <- list(
        Sigma_raw = Sigma_raw_grp,
        Sigma_multilag = Sigma_ml_grp,
        R_group = R_grp,
        C_sym = C_u,
        M_op = sm_grp$M_op,
        eta_opt = sm_grp$eta_opt
      )
    }
  }

  structure(
    list(
      coefficients = data.table::rbindlist(coef_out_list, fill = TRUE),
      latency = data.table::rbindlist(lat_out_list, fill = TRUE),
      contrasts = if (length(con_out_list)) data.table::rbindlist(con_out_list, fill = TRUE) else data.table::data.table(),
      subjects = if (length(subj_out_list)) data.table::rbindlist(subj_out_list, fill = TRUE) else NULL,
      subject_latency = if (length(subj_lat_list)) data.table::rbindlist(subj_lat_list, fill = TRUE) else NULL,
      covariance = cov_out_list,
      settings = list(
        lags = lags,
        trajectory_mode = trajectory_mode,
        demix_within_trial = demix_within_trial,
        n_perm = n_perm,
        n_boot_latency = n_boot_latency,
        fine_dt = fine_dt,
        conf_level = conf_level,
        polarity = polarity
      )
    ),
    class = "medusa_multilag"
  )
}


