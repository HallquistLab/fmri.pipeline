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
