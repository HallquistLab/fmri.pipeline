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
  quad <- epoch_kernel_quadrature(readout_kernel, readout_width, tr)
  m <- length(onsets)
  readout <- matrix(0, m, n_time)
  weight_sum <- rep(0, m)
  for (q in seq_along(quad$offsets)) {
    at <- (onsets + lag + quad$offsets[q] - offset) / tr
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
#'   \code{"impulse"} (default unit impulse at \code{onset + lag}) or
#'   \code{"window"} (local sustained boxcar of width \code{readout_width}
#'   centered at \code{onset + lag}).
#' @param return_epochs Logical; whether to include raw and de-mixed trial-by-voxel
#'   epoch matrices per lag. Default: \code{FALSE}.
#' @param return_panel Logical; whether to return a tidy MEDUSA-compatible
#'   \code{data.table} (\code{panel}) with columns \code{run}, \code{trial},
#'   \code{event_index}, \code{onset}, \code{evt_time}, \code{voxel},
#'   \code{decon_raw}, \code{decon_demixed}, and the covariate columns for
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
  source_model = c("impulse", "window"),
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
  readout_kernel <- match.arg(readout_kernel)
  source_model <- match.arg(source_model)
  checkmate::assert_numeric(lags, min.len = 1L, any.missing = FALSE, finite = TRUE)
  checkmate::assert_string(target_condition, null.ok = TRUE)
  checkmate::assert_flag(demix_other_events)
  checkmate::assert_flag(demix_target_trials)
  checkmate::assert_number(demix_ridge, lower = 1e-6, finite = TRUE)
  checkmate::assert_number(gain_floor, lower = 1e-4, finite = TRUE)
  checkmate::assert_flag(restore_scale)
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

  raw_epochs <- demixed_epochs <- vector("list", n_lags)
  names(raw_epochs) <- names(demixed_epochs) <- sprintf("lag_%g", lags)
  for (li in seq_len(n_lags)) {
    raw_epochs[[li]] <- matrix(NA_real_, sum(target_counts), n_voxels,
      dimnames = list(NULL, fit$voxel_ids))
    demixed_epochs[[li]] <- matrix(NA_real_, sum(target_counts), n_voxels,
      dimnames = list(NULL, fit$voxel_ids))
  }
  gain_by_lag_voxel <- matrix(NA_real_, n_lags, n_voxels, dimnames = list(NULL, fit$voxel_ids))
  x_centered_list <- lapply(cov_list, function(x_r) sweep(x_r, 2L, colMeans(x_r, na.rm = TRUE), "-"))

  ends <- cumsum(target_counts)
  starts <- ends - target_counts + 1L
  gain_accum <- matrix(0, n_lags, n_voxels)
  gain_count <- matrix(0L, n_lags, n_voxels)

  for (r in seq_len(n_runs)) {
    n <- fit$n_time[r]
    ev <- events[[r]]
    offset <- fit$settings$time_offset[r]
    t_idx <- target_idx_list[[r]]
    o_idx <- other_idx_list[[r]]
    m_all <- nrow(ev)
    m_target <- length(t_idx)
    row_span <- starts[r]:ends[r]
    duration <- if (is.null(ev$duration)) rep(0, m_all) else ev$duration

    mask <- if (settings$observation_mask == "fit") fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    if (length(rows) < 3L) next
    z <- if (settings$nuisance_mode == "joint") {
      cbind(1, seq(-1, 1, length.out = n), fit$specification$nuisance[[r]])
    } else matrix(numeric(), n, 0L)

    selected <- fit$selections[[r]]$kernel_index
    rho <- if (is.null(fit$noise)) rep(0, length(selected)) else fit$noise[[r]]$rho
    keys <- paste(selected, sprintf("%.17g", rho), sep = ":")
    valid <- deconvolution$scaling[[r]]$status == "ok"
    act_r <- deconvolution$activity[[r]]
    if (isTRUE(restore_scale)) {
      div_r <- deconvolution$scaling[[r]]$scale
      div_r[!is.finite(div_r) | div_r <= 0] <- 1
      act_r[, valid] <- sweep(act_r[, valid, drop = FALSE], 2L, div_r[valid], "*")
    }

    S_all_0 <- build_epoch_readout_matrix(ev$onset, 0, n, tr, offset,
      readout_kernel = readout_kernel, readout_width = readout_width)
    S_lags <- lapply(lags, function(lag_val) {
      build_epoch_readout_matrix(ev$onset[t_idx], lag_val, n, tr, offset,
        readout_kernel = readout_kernel, readout_width = readout_width)
    })

    for (li in seq_len(n_lags)) {
      raw_epochs[[li]][row_span, valid] <- S_lags[[li]] %*% act_r[, valid, drop = FALSE]
      demixed_epochs[[li]][row_span, valid] <- raw_epochs[[li]][row_span, valid, drop = FALSE]
    }

    for (key in unique(keys[valid])) {
      voxels <- which(keys == key & valid)
      k <- selected[voxels[1L]]
      ar <- rho[voxels[1L]]
      hrf_k <- fit$kernels[, k, drop = FALSE]
      operator <- hrf_inverse_operator(n, hrf_k[, 1L], rows, ar, z, settings)

      # Forward design for all events at their own onsets (for non-target leakage removal)
      G_all_0 <- hrf_event_design_cpp(n, tr, offset, ev$onset, duration,
        seq_len(m_all), m_all, hrf_k)[[1L]][rows, , drop = FALSE]
      if (settings$scaling != "none") G_all_0 <- sweep(G_all_0, 2L, colMeans(G_all_0), "-")
      DG_all_0 <- operator %*% G_all_0

      a_other_hat <- NULL
      if (demix_other_events && length(o_idx) > 0L) {
        ok_all <- rowSums(!is.finite(S_all_0)) == 0L
        if (any(ok_all)) {
          # Estimate condition-level mean amplitudes of all event types jointly, plus
          # stabilized trial-specific residuals of non-target events orthogonal to target events
          cond_all <- if (is.null(ev$condition)) rep("event", m_all) else as.character(ev$condition)
          cond_levels <- unique(cond_all[ok_all])
          M_cond <- vapply(cond_levels, function(cl) as.numeric(cond_all == cl), numeric(m_all))
          if (is.null(dim(M_cond))) M_cond <- matrix(M_cond, ncol = 1L)
          R_all_0 <- S_all_0[ok_all, , drop = FALSE] %*% DG_all_0
          Y_all_0 <- S_all_0[ok_all, , drop = FALSE] %*% act_r[, voxels, drop = FALSE]
          R_cond <- R_all_0 %*% M_cond
          lhs_cond <- crossprod(R_cond) + 1e-4 * mean(diag(crossprod(R_cond))) * diag(ncol(R_cond))
          alpha_cond <- solve(lhs_cond, crossprod(R_cond, Y_all_0))
          a_mean_full <- M_cond %*% alpha_cond

          # Residual at non-target onsets after removing condition means, projected
          # orthogonal to the target events' resolution subspace so target signal is never subtracted
          ok_o <- o_idx[o_idx %in% which(ok_all)]
          if (length(ok_o)) {
            R_oo <- S_all_0[ok_o, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
            R_ot <- S_all_0[ok_o, , drop = FALSE] %*% DG_all_0[, t_idx, drop = FALSE]
            Y_res_o <- (S_all_0[ok_o, , drop = FALSE] %*% act_r[, voxels, drop = FALSE]) -
              (S_all_0[ok_o, , drop = FALSE] %*% DG_all_0) %*% a_mean_full
            if (ncol(R_ot) > 0L) {
              diag_ot <- max(mean(diag(crossprod(R_ot))), 1e-6)
              proj_t <- solve(crossprod(R_ot) + 0.05 * diag_ot * diag(ncol(R_ot)), crossprod(R_ot, Y_res_o))
              Y_res_o <- Y_res_o - R_ot %*% proj_t
            }
            diag_oo <- max(mean(diag(crossprod(R_oo))), 1e-6)
            lhs_oo <- crossprod(R_oo) + max(demix_ridge, 0.20) * diag_oo * diag(ncol(R_oo))
            delta_o <- solve(lhs_oo, crossprod(R_oo, Y_res_o))
            a_other_hat <- a_mean_full[o_idx, , drop = FALSE] + delta_o
          } else {
            a_other_hat <- a_mean_full[o_idx, , drop = FALSE]
          }
        }
      }

      for (li in seq_len(n_lags)) {
        S_tau <- S_lags[[li]]
        ok_tau <- rowSums(!is.finite(S_tau)) == 0L
        if (!any(ok_tau)) next
        Y_clean <- raw_epochs[[li]][row_span, voxels, drop = FALSE]

        if (!is.null(a_other_hat)) {
          R_other_tau <- S_tau[ok_tau, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
          Y_clean[ok_tau, ] <- Y_clean[ok_tau, , drop = FALSE] - R_other_tau %*% a_other_hat
        }

        # Local target resolution matrix K_tau for unit sources at target onset + lag
        w_src <- if (is.null(readout_width)) 1.5 * tr else readout_width
        dur_tau <- if (source_model == "window") rep(w_src, m_target) else duration[t_idx]
        shift_center <- if (source_model == "window") lags[li] - 0.5 * w_src else lags[li]
        shifted_onsets <- pmin(max(ev$onset), pmax(offset, ev$onset[t_idx] + shift_center))
        G_tau <- hrf_event_design_cpp(n, tr, offset, shifted_onsets, dur_tau,
          seq_len(m_target), m_target, hrf_k)[[1L]][rows, , drop = FALSE]
        if (settings$scaling != "none") G_tau <- sweep(G_tau, 2L, colMeans(G_tau), "-")
        K_tau <- S_tau[ok_tau, , drop = FALSE] %*% (operator %*% G_tau[, ok_tau, drop = FALSE])

        # For positive-duration boxcar events or window source_model, normalize by ideal readout diagonal;
        # for impulse events (duration == 0), G_tau is already a unit impulse source.
        if (any(dur_tau > 0)) {
          neural_tau <- hrf_event_design_cpp(n, tr, offset, shifted_onsets, dur_tau,
            seq_len(m_target), m_target, matrix(c(1, 0), ncol = 1L))[[1L]]
          ideal_diag <- diag(S_tau[ok_tau, , drop = FALSE] %*% neural_tau[, ok_tau, drop = FALSE])
          ideal_diag[!is.finite(ideal_diag) | abs(ideal_diag) < 1e-4] <- 1
          d_pos <- dur_tau[ok_tau] > 0
          K_norm <- K_tau
          if (any(d_pos)) K_norm[d_pos, ] <- K_tau[d_pos, , drop = FALSE] / ideal_diag[d_pos]
        } else {
          K_norm <- K_tau
        }
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
          Y_demix_unit <- solve(lhs_c, crossprod(C_tau, Y_clean[ok_tau, , drop = FALSE]))
          Y_clean[ok_tau, ] <- Y_demix_unit / g_safe
        }
        demixed_epochs[[li]][row_span, voxels] <- Y_clean
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
      ok_v <- is.finite(y_raw_v) & is.finite(y_dem_v) & (rowSums(!is.finite(X_pool)) == 0L)
      n_ok <- sum(ok_v)
      if (n_ok <= n_terms + 2L) {
        for (ti in seq_len(n_terms)) {
          row_k <- row_k + 1L
          coef_rows[[row_k]] <- data.frame(
            voxel = fit$voxel_ids[vi], lag = lags[li], term = term_names[ti],
            beta_raw = NA_real_, beta_corrected = NA_real_,
            r_raw = NA_real_, r_demixed = NA_real_,
            t_stat = NA_real_, p_value = NA_real_, df = n_ok - n_terms,
            mean_gain = gain_by_lag_voxel[li, vi], stringsAsFactors = FALSE
          )
        }
        next
      }

      X_ok <- X_pool[ok_v, , drop = FALSE]
      xtx_raw <- crossprod(X_ok) + 1e-8 * diag(n_terms)
      beta_raw_v <- as.numeric(solve(xtx_raw, crossprod(X_ok, y_raw_v[ok_v])))
      beta_dem_v <- as.numeric(solve(xtx_raw, crossprod(X_ok, y_dem_v[ok_v])))

      resid_dem <- y_dem_v[ok_v] - as.numeric(X_ok %*% beta_dem_v)
      df_v <- max(1L, n_ok - n_terms - n_runs)
      sigma2_dem <- sum(resid_dem^2) / df_v
      se_dem <- sqrt(pmax(diag(solve(xtx_raw)) * sigma2_dem, .Machine$double.eps))
      t_vec <- beta_dem_v / se_dem
      p_vec <- 2 * stats::pt(abs(t_vec), df = df_v, lower.tail = FALSE)

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
          beta_raw = beta_raw_v[ti], beta_corrected = beta_dem_v[ti],
          r_raw = r_raw_val, r_demixed = r_dem_val,
          t_stat = t_vec[ti], p_value = p_vec[ti], df = df_v,
          mean_gain = gain_by_lag_voxel[li, vi], stringsAsFactors = FALSE
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
        mean_r_raw = fmean(sub_df$r_raw),
        mean_r_demixed = fmean(sub_df$r_demixed),
        mean_gain = fmean(sub_df$mean_gain),
        n_voxels = sum(is.finite(sub_df$beta_corrected)),
        stringsAsFactors = FALSE
      )
    }
  }

  panel_dt <- NULL
  if (isTRUE(return_panel)) {
    total_trials <- sum(target_counts)
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
        decon_demixed = as.numeric(demixed_epochs[[li]])
      )
      for (cn in colnames(cov_full)) {
        base_meta[[cn]] <- rep(cov_full[[cn]], times = n_voxels)
      }
      panel_list[[li]] <- base_meta
    }
    panel_dt <- data.table::rbindlist(panel_list)
  }

  list(
    coefficients = coef_df,
    summary = do.call(rbind, summary_rows),
    epochs = if (return_epochs) list(raw = raw_epochs, demixed = demixed_epochs) else NULL,
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
      source_model = source_model
    )
  )
}
