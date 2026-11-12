#' Diagnose temporal spread of a fixed linear deconvolution operator
#'
#' Generate isolated unit impulses at actual event times, convolve them with a
#' known HRF, and measure their recovery and entry into other event windows.
#'
#' @param event_data Data.frame with unique trial labels and numeric onsets.
#' @param kernel Frozen inverse-model HRF, normalized according to
#'   \code{kernel_normalize}.
#' @param n_time Number of measured acquisition samples (at least three).
#' @param tr Acquisition interval in seconds.
#' @param lambda Fixed regularization strength; automatic tuning is not performed.
#' @param onset,trial Names of onset and trial columns in event_data.
#' @param time_offset Absolute time assigned to the first acquisition.
#' @param time_before,time_after Requested window limits relative to each event.
#' @param output_resolution Spacing of the diagnostic window grid, in seconds.
#' @param backend Linear backend: cpp, r, or scalar wiener. Adaptive PSD and
#'   nonlinear estimators are not supported.
#' @param penalty Fixed quadratic penalty for cpp/r.
#' @param forward_kernel Generating HRF, defaulting to kernel. Use a different
#'   shape for a specified HRF-mismatch experiment, normalized with the same
#'   \code{kernel_normalize} convention as \code{kernel}.
#' @param prewhiten_rho Frozen AR(1) coefficient; zero disables whitening.
#' @param ridge_floor Fixed numerical ridge floor for cpp/r.
#' @param input_scale Fixed positive BOLD divisor from the fitted run. Recovery is
#'   multiplied back by this divisor to report activity in injected-impulse units.
#' @param input_demean Whether to subtract each simulated measured BOLD mean
#'   before padding. This applies a fixed linear centering operation, not z-scoring.
#' @param tail_padding Number of zero BOLD samples appended as observations.
#'   Use length(kernel) to reproduce the current voxelwise wrapper's padding.
#' @param source_trials Optional trial labels to inject; all are used by default.
#'   All event_data rows remain target windows. Subset sources for large schedules.
#' @param batch_size Maximum number of impulse columns fitted together.
#' @param return_timeseries Include full acquisition-grid truth, BOLD, and
#'   recovered activity for each selected source (default FALSE).
#' @param kernel_normalize Normalization applied to both \code{kernel} and
#'   \code{forward_kernel}: \code{"auto"} (default; honors \code{attr(..., "normalize")}
#'   when present, switches to \code{"positive_lobe"} when \code{kernel} and
#'   \code{forward_kernel} differ in undershoot ratio, and uses
#'   \code{"signed_sum"} otherwise), \code{"signed_sum"}, \code{"positive_lobe"},
#'   or \code{"peak"}.
#' @return A list with impulses (one status/summary per selected source), windows
#'   (directed source/target metrics), optional time_series, and frozen settings.
#' @details Off-grid events deposit mass on the two surrounding acquisition bins;
#'   this is a grid approximation, not a sub-TR neural-resolution claim. The
#'   injected samples sum to one (area tr). Missing/outside onsets and impulses
#'   producing no observed BOLD variation retain explicit statuses and missing
#'   recovery metrics. Windows use overlapping support without collision censoring.
#'   window_l1_ratio compares absolute recovered area in a target window with that
#'   in the source's own window. Windows may overlap, so ratios are not fractions
#'   of neural attribution and may exceed one. Metrics condition on fixed model
#'   settings and the specified generating HRF. No noise, retuning, observation-
#'   mask refitting, or physiological uncertainty estimation is performed.
#' @export
deconvolution_impulse_diagnostics <- function(event_data, kernel, n_time, tr, lambda,
    onset = "onset", trial = "trial", time_offset = 0, time_before = -3, time_after = 3,
    output_resolution = tr, backend = c("cpp", "r", "wiener"),
    penalty = c("diff2", "diff1", "ridge", "matern1", "sobolev"), forward_kernel = NULL,
    prewhiten_rho = 0, ridge_floor = 1e-8, input_scale = 1, input_demean = FALSE,
    tail_padding = 0L, source_trials = NULL, batch_size = 32L, return_timeseries = FALSE,
    kernel_normalize = c("auto", "signed_sum", "positive_lobe", "peak")) {
  backend <- match.arg(backend)
  penalty <- match.arg(penalty)
  kernel_normalize <- match.arg(kernel_normalize)
  checkmate::assert_data_frame(event_data, min.rows = 1L)
  checkmate::assert_string(onset)
  checkmate::assert_string(trial)
  checkmate::assert_subset(c(onset, trial), names(event_data))
  checkmate::assert_numeric(event_data[[onset]], all.missing = TRUE)
  labels <- as.character(event_data[[trial]])
  if (anyNA(labels) || anyDuplicated(labels) || any(!nzchar(labels))) stop("Trial labels must be unique and nonmissing")
  checkmate::assert_integerish(n_time, lower = 3L, len = 1L)
  checkmate::assert_number(tr, lower = 0.01, finite = TRUE)
  checkmate::assert_number(lambda, lower = 0, finite = TRUE)
  checkmate::assert_number(time_offset, finite = TRUE)
  checkmate::assert_number(prewhiten_rho, lower = -0.99, upper = 0.99, finite = TRUE)
  checkmate::assert_number(ridge_floor, lower = 0, finite = TRUE)
  checkmate::assert_number(input_scale, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_flag(input_demean)
  checkmate::assert_integerish(tail_padding, lower = 0, len = 1L)
  checkmate::assert_integerish(batch_size, lower = 1L, len = 1L)
  checkmate::assert_flag(return_timeseries)
  checkmate::assert_numeric(kernel, min.len = 2L, any.missing = FALSE, finite = TRUE)
  if (!is.null(forward_kernel)) {
    checkmate::assert_numeric(forward_kernel, min.len = 2L, any.missing = FALSE, finite = TRUE)
  }
  norm_mode <- if (kernel_normalize != "auto") {
    kernel_normalize
  } else if (!is.null(attr(kernel, "normalize"))) {
    match.arg(attr(kernel, "normalize"), c("signed_sum", "positive_lobe", "peak"))
  } else if (!is.null(forward_kernel) && !is.null(attr(forward_kernel, "normalize"))) {
    match.arg(attr(forward_kernel, "normalize"), c("signed_sum", "positive_lobe", "peak"))
  } else if (!is.null(forward_kernel)) {
    infer_kernel_grid_normalize(list(kernel, forward_kernel))
  } else {
    "signed_sum"
  }
  kernel <- normalize_decon_kernel(kernel, mode = norm_mode)
  forward_kernel <- if (is.null(forward_kernel)) kernel else normalize_decon_kernel(forward_kernel, mode = norm_mode)
  grid <- epoch_output_grid(time_before, time_after, output_resolution)
  source <- if (is.null(source_trials)) seq_along(labels) else {
    source_trials <- as.character(source_trials)
    checkmate::assert_character(source_trials, min.len = 1L, unique = TRUE, any.missing = FALSE)
    checkmate::assert_subset(source_trials, labels)
    match(source_trials, labels)
  }
  n_time <- as.integer(n_time)
  time <- time_offset + (seq_len(n_time) - 1L) * tr
  onsets <- event_data[[onset]]
  position <- (onsets - time_offset) / tr
  position[is.finite(position) & abs(position) < 1e-8] <- 0
  position[is.finite(position) & abs(position - (n_time - 1L)) < 1e-8] <- n_time - 1L
  target_times <- outer(onsets, grid, "+")
  status <- ifelse(!is.finite(position), "missing_onset",
    ifelse(position < 0 | position > n_time - 1L, "outside_run", "ok"))
  # This causal forward operator has zero pre-run neural history.
  forward <- decon_convolution_matrix(n_time, forward_kernel)
  summaries <- windows <- series <- list()
  batches <- split(source, ceiling(seq_along(source) / batch_size))
  for (indices in batches) {
    truth <- matrix(0, n_time, length(indices))
    for (j in seq_along(indices)) {
      i <- indices[j]
      if (status[i] != "ok") next
      lo <- floor(position[i])
      weight <- position[i] - lo
      truth[lo + 1L, j] <- 1 - weight
      if (weight > 0) truth[lo + 2L, j] <- weight
    }
    bold <- forward %*% rbind(matrix(0, length(forward_kernel) - 1L, ncol(truth)), truth)
    y <- if (input_demean) sweep(bold, 2L, colMeans(bold), "-") else bold
    y <- y / input_scale
    valid <- status[indices] == "ok" & decon_signal_status(y) == "ok"
    recovered <- matrix(NA_real_, n_time, length(indices))
    if (any(valid)) {
      fit_input <- rbind(y[, valid, drop = FALSE], matrix(0, tail_padding, sum(valid)))
      recovered[, valid] <- deconvolve_reglin(fit_input, kernel, lambda = lambda,
        backend = backend, penalty = penalty, scaling = "none", ridge_floor = ridge_floor,
        prewhiten_gcv = prewhiten_rho != 0, prewhiten_rho = prewhiten_rho,
        kernel_normalize = norm_mode)[seq_len(n_time), , drop = FALSE] * input_scale
    }
    for (j in seq_along(indices)) {
      i <- indices[j]
      source_status <- status[i]
      if (source_status == "ok" && !valid[j]) {
        source_status <- if (all(bold[, j] == 0)) "unobserved_impulse" else "constant_bold"
      }
      activity <- recovered[, j]
      estimate <- function(at) {
        interpolate_epoch_unit(time, activity, seq_len(n_time), rep(FALSE, n_time), at, tr)$.estimate
      }
      target_values <- matrix(NA_real_, nrow(target_times), ncol(target_times))
      finite_targets <- is.finite(target_times)
      if (any(is.finite(activity)) && any(finite_targets)) {
        target_values[finite_targets] <- estimate(target_times[finite_targets])
      }
      source_windows <- lapply(seq_along(labels), function(target) {
        values <- target_values[target, ]
        good <- is.finite(values)
        data.frame(source_trial = labels[i], target_trial = labels[target],
          source_onset = onsets[i], target_onset = onsets[target], source_status = source_status,
          target_onset_status = status[target], n_lags = length(grid), n_valid = sum(good),
          window_l1 = if (any(good)) sum(abs(values[good])) * output_resolution else NA_real_,
          window_signed_area = if (any(good)) sum(values[good]) * output_resolution else NA_real_,
          window_max_abs = if (any(good)) max(abs(values[good])) else NA_real_)
      })
      source_windows <- do.call(rbind, source_windows)
      denominator <- source_windows$window_l1[i]
      source_windows$window_l1_ratio <- if (is.finite(denominator) && denominator > 0) {
        source_windows$window_l1 / denominator
      } else NA_real_
      windows[[length(windows) + 1L]] <- source_windows
      available <- any(is.finite(activity))
      peak <- if (available) which.max(abs(activity)) else NA_integer_
      truth_onset <- if (status[i] == "ok") interpolate_epoch_unit(time, truth[, j],
        seq_len(n_time), rep(FALSE, n_time), onsets[i], tr)$.estimate else NA_real_
      at_onset <- if (status[i] == "ok") estimate(onsets[i]) else NA_real_
      summaries[[length(summaries) + 1L]] <- data.frame(source_trial = labels[i],
        source_onset = onsets[i], status = source_status, truth_at_onset = truth_onset,
        estimate_at_onset = at_onset, onset_gain = at_onset / truth_onset,
        peak_amplitude = if (available) activity[peak] else NA_real_,
        peak_time_relative = if (available) time[peak] - onsets[i] else NA_real_,
        mass_gain = if (available) sum(activity) else NA_real_,
        recovered_l1 = if (available) sum(abs(activity)) * tr else NA_real_,
        pre_onset_l1 = if (available) sum(abs(activity[time < onsets[i]])) * tr else NA_real_,
        negative_mass = if (available) sum(pmax(-activity, 0)) * tr else NA_real_,
        distance_run_start = onsets[i] - min(time), distance_run_end = max(time) - onsets[i])
      if (return_timeseries) {
        series[[length(series) + 1L]] <- data.frame(source_trial = labels[i], time = time,
          relative_time = time - onsets[i], truth = if (status[i] == "ok") truth[, j] else NA_real_,
          bold = if (status[i] == "ok") bold[, j] else NA_real_, activity = activity)
      }
    }
  }
  list(impulses = do.call(rbind, summaries), windows = do.call(rbind, windows),
    time_series = if (return_timeseries) do.call(rbind, series) else NULL,
    settings = list(kernel = kernel, forward_kernel = forward_kernel,
      kernel_normalize = norm_mode, lambda = lambda,
      backend = backend, penalty = penalty, n_time = n_time, tr = tr, time_offset = time_offset,
      time_before = time_before, time_after = time_after, output_resolution = output_resolution,
      prewhiten_rho = prewhiten_rho, ridge_floor = ridge_floor, input_scale = input_scale,
      input_demean = input_demean, tail_padding = tail_padding,
      onset_deposition = "linear_acquisition_bins", attribution = FALSE))
}
