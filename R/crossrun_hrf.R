#' Build a small, constrained library of causal HRFs
#'
#' @param tr Acquisition interval in seconds.
#' @param peak_delays Positive response-delay parameters in seconds.
#' @param dispersions Positive response-dispersion parameters.
#' @param undershoot_ratios Positive response-to-undershoot ratios.
#' @param duration Kernel duration in seconds (at least two acquisition intervals).
#' @return A unit-sum kernel matrix, with time rows, candidate columns, and a
#'   parameters attribute. The canonical candidate (6, 1, 6) is first if present.
#' @details Candidates use the SPM difference-of-gammas parameterization with
#'   undershoot delay 16 seconds and undershoot dispersion 1. These are candidate
#'   shapes, not a universally validated range for every population or protocol.
#' @export
hrf_library <- function(tr, peak_delays = c(5, 6, 7), dispersions = c(0.8, 1, 1.2),
    undershoot_ratios = 6, duration = 32) {
  checkmate::assert_number(tr, lower = 0.01, finite = TRUE)
  for (x in list(peak_delays, dispersions, undershoot_ratios))
    checkmate::assert_numeric(x, lower = .Machine$double.eps, min.len = 1L,
      any.missing = FALSE, finite = TRUE, unique = TRUE)
  checkmate::assert_number(duration, lower = 2 * tr, finite = TRUE)
  pars <- expand.grid(peak_delay = peak_delays, dispersion = dispersions,
    undershoot_ratio = undershoot_ratios, KEEP.OUT.ATTRS = FALSE)
  canonical <- with(pars, peak_delay == 6 & dispersion == 1 & undershoot_ratio == 6)
  pars <- pars[order(!canonical, seq_len(nrow(pars))), , drop = FALSE]
  rownames(pars) <- NULL
  kernels <- vapply(seq_len(nrow(pars)), function(i) {
    as.numeric(spm_hrf(tr, c(pars$peak_delay[i], 16, pars$dispersion[i], 1,
      pars$undershoot_ratio[i], 0, duration))$hrf)
  }, numeric(floor(duration / tr) + 1L))
  colnames(kernels) <- paste0("hrf_", seq_len(ncol(kernels)))
  attr(kernels, "parameters") <- pars
  attr(kernels, "tr") <- tr
  kernels
}

hrf_validate_bold <- function(bold) {
  checkmate::assert_list(bold, min.len = 2L, names = "unique")
  if (any(!nzchar(names(bold)))) stop("BOLD runs must have nonempty unique names")
  for (y in bold) {
    checkmate::assert_matrix(y, mode = "numeric", min.rows = 3L, min.cols = 1L)
    if (is.null(colnames(y)) || anyNA(colnames(y)) || any(!nzchar(colnames(y))) ||
        anyDuplicated(colnames(y))) stop("BOLD matrices require unique voxel column names")
    if (!identical(colnames(y), colnames(bold[[1L]])))
      stop("All BOLD runs must have identical voxel IDs in the same order")
  }
  invisible(TRUE)
}

hrf_run_list <- function(x, bold, argument, default) {
  if (is.null(x)) return(lapply(bold, default))
  if (!is.list(x) || !identical(names(x), names(bold)))
    stop(argument, " must be a list with the same run names and order as bold")
  x
}

hrf_run_numeric <- function(x, bold, argument) {
  checkmate::assert_numeric(x, any.missing = FALSE, finite = TRUE)
  if (length(x) == 1L) return(stats::setNames(rep(x, length(bold)), names(bold)))
  if (length(x) != length(bold) || !identical(names(x), names(bold)))
    stop(argument, " must be scalar or named in BOLD run order")
  x
}

#' Estimate HRF shapes from other runs, with optional spatial regularization
#'
#' @param bold Named list of time-by-voxel numeric matrices for one participant.
#'   Voxel column IDs must be unique and identical in all runs, in the same order.
#' @param events Named list of data.frames with onset (seconds), condition, and
#'   optional duration (seconds, default zero). Names/order must match bold.
#' @param tr Common acquisition interval in seconds. All runs use the same grid.
#' @param kernel_grid Matrix or list of equal-length candidate causal HRFs, sampled
#'   at tr starting at lag zero. Use hrf_library() or supply explicit candidates.
#' @param nuisance Optional matching list of time-by-regressor matrices. An
#'   intercept and linear trend are always included separately for each run.
#' @param valid_rows Optional matching list of logical observation masks; FALSE
#'   and NA omit rows from training fits without deleting their acquisition times.
#' @param time_offset First acquisition time, scalar or named numeric run vector.
#' @param spatial_groups Optional voxel parcel labels, named in voxel column order.
#'   Labels must be nonmissing and nonempty. Parcels must be defined independently
#'   of target-run responses; no implicit whole-brain parcel is constructed.
#' @param pooling Select HRFs per voxel or share one HRF per supplied parcel.
#' @param spatial_strength Nonnegative fixed penalty on voxel candidate distance
#'   from the training-derived parcel HRF; requires spatial_groups. Zero disables
#'   the soft penalty. Only used with voxel pooling; not automatically tuned.
#' @param fallback_kernel Index of the prespecified fallback candidate, default 1.
#' @param min_r2 Minimum mean adjusted partial R-squared for a supported selection,
#'   default 0.01. This is a configurable heuristic, not a significance threshold.
#' @param min_training_runs Minimum usable other runs per voxel, default 1.
#' @param min_parcel_voxels Minimum eligible training voxels per parcel, default 5.
#' @param batch_size Maximum voxel columns per native scoring block.
#' @param return_scores Retain candidate-by-voxel training score matrices per fold.
#' @param amplitude_model Training amplitudes shared within condition (fast
#'   default), or free for each event record (trial). Trial amplitudes can support
#'   HRF estimation when signed trial responses average to zero, but cost more and
#'   require a design with sufficient residual degrees of freedom.
#' @return An object of class crossrun_hrf with normalized kernels, per-target-run
#'   selections and training manifests, parcel anchors, run status/design summaries,
#'   optional candidate scores, voxel/run IDs, and frozen settings.
#' @details Each target run is excluded from its own HRF selection, including
#'   parcel anchors and fallback decisions. Training fits have run-specific
#'   condition or event amplitudes and nuisance coefficients, with one shared
#'   candidate HRF across conditions. They do not fit unrestricted trial
#'   timecourses or use RPE.
#'   Native SVD projections are shared across voxel blocks. Scores are adjusted
#'   partial R-squared after removing nuisance, averaged equally across usable
#'   training runs. Parcel scores average these voxel scores equally, retaining
#'   voxel-specific amplitudes rather than averaging BOLD first. In voxel mode,
#'   selection maximizes score minus spatial_strength times normalized squared
#'   distance to the reliable parcel anchor. Weak/insufficient voxel fits fall
#'   back to a reliable parcel anchor, then the prespecified fallback kernel.
#'   This is a discrete shape penalty, not Bayesian posterior inference, a local
#'   adjacency prior, or spatial smoothing of neural activity. Library, thresholds,
#'   parcels, and penalty must be fixed externally or tuned within training data.
#'   No automatic nested tuning, uncertainty intervals, or AR noise fit is included.
#'   Nonfinite/constant/nuisance-only voxel runs are omitted; observation masks are
#'   common across voxels. Candidate task ranks must match within a run and leave
#'   residual degrees of freedom, otherwise that run cannot inform HRF selection.
#' @export
estimate_crossrun_hrf <- function(bold, events, tr, kernel_grid, nuisance = NULL,
    valid_rows = NULL, time_offset = 0, spatial_groups = NULL,
    pooling = c("voxel", "parcel"), spatial_strength = 0, fallback_kernel = 1L,
    min_r2 = 0.01, min_training_runs = 1L, min_parcel_voxels = 5L,
    batch_size = 2048L, return_scores = FALSE, amplitude_model = c("condition", "trial")) {
  hrf_validate_bold(bold)
  pooling <- match.arg(pooling)
  amplitude_model <- match.arg(amplitude_model)
  checkmate::assert_number(tr, lower = 0.01, finite = TRUE)
  if (!is.null(attr(kernel_grid, "tr")) && !isTRUE(all.equal(attr(kernel_grid, "tr"), tr)))
    stop("kernel_grid sampling TR does not match tr")
  parameters <- attr(kernel_grid, "parameters")
  kernels <- build_decon_kernel_grid(c(0, 1), kernel_grid)
  if (any(vapply(kernels, function(h) sum(abs(h)) > 100, logical(1))))
    stop("Candidate signed sums are too close to zero for stable unit-sum normalization")
  kernels <- do.call(cbind, kernels)
  candidate_ids <- if (is.matrix(kernel_grid)) colnames(kernel_grid) else names(kernel_grid)
  if (is.null(candidate_ids)) candidate_ids <- paste0("hrf_", seq_len(ncol(kernels)))
  if (anyNA(candidate_ids) || any(!nzchar(candidate_ids)) || anyDuplicated(candidate_ids))
    stop("Candidate names must be unique and nonempty")
  colnames(kernels) <- candidate_ids
  checkmate::assert_integerish(fallback_kernel, lower = 1, upper = ncol(kernels), len = 1L)
  checkmate::assert_number(spatial_strength, lower = 0, finite = TRUE)
  checkmate::assert_number(min_r2, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_integerish(min_training_runs, lower = 1, upper = length(bold) - 1L, len = 1L)
  checkmate::assert_integerish(min_parcel_voxels, lower = 1, len = 1L)
  checkmate::assert_integerish(batch_size, lower = 1, len = 1L)
  checkmate::assert_flag(return_scores)
  voxel_ids <- colnames(bold[[1L]])
  if (!is.null(spatial_groups)) {
    if (!is.atomic(spatial_groups) || !identical(names(spatial_groups), voxel_ids) ||
        anyNA(spatial_groups) || any(!nzchar(as.character(spatial_groups))))
      stop("spatial_groups must be named nonmissing labels in voxel column order")
    spatial_groups <- as.character(spatial_groups)
  } else if (pooling == "parcel" || spatial_strength > 0) {
    stop("Parcel pooling/spatial regularization requires explicit spatial_groups")
  }
  if (pooling == "parcel" && spatial_strength != 0)
    stop("spatial_strength applies only to voxel pooling")
  offsets <- hrf_run_numeric(time_offset, bold, "time_offset")
  events <- hrf_run_list(events, bold, "events", function(y) NULL)
  nuisance <- hrf_run_list(nuisance, bold, "nuisance", function(y) matrix(numeric(), nrow(y), 0))
  valid_rows <- hrf_run_list(valid_rows, bold, "valid_rows", function(y) rep(TRUE, nrow(y)))
  per_run <- design_manifest <- event_manifest <- vector("list", length(bold))
  names(per_run) <- names(design_manifest) <- names(event_manifest) <- names(bold)
  if (!is.null(spatial_groups)) names(spatial_groups) <- voxel_ids
  for (r in seq_along(bold)) {
    y <- bold[[r]]; ev <- events[[r]]
    checkmate::assert_data_frame(ev, min.rows = 1L)
    checkmate::assert_subset(c("onset", "condition"), names(ev))
    checkmate::assert_numeric(ev$onset, any.missing = FALSE, finite = TRUE)
    labels <- as.character(ev$condition)
    if (anyNA(labels) || any(!nzchar(labels))) stop("Event conditions must be nonmissing and nonempty")
    durations <- if (is.null(ev$duration)) rep(0, nrow(ev)) else ev$duration
    checkmate::assert_numeric(durations, lower = 0, any.missing = FALSE, finite = TRUE)
    checkmate::assert_matrix(nuisance[[r]], mode = "numeric", nrows = nrow(y))
    mask <- valid_rows[[r]]
    checkmate::assert_logical(mask, len = nrow(y))
    mask <- !is.na(mask) & mask
    design_columns <- if (amplitude_model == "trial") seq_len(nrow(ev)) else match(labels, unique(labels))
    design <- hrf_event_design_cpp(nrow(y), tr, offsets[r], ev$onset, durations,
      design_columns, max(design_columns), kernels)
    z <- cbind(intercept = 1, trend = seq(-1, 1, length.out = nrow(y)), nuisance[[r]])
    result <- if (sum(mask) >= 3L) {
      hrf_score_run_cpp(y, design, z, which(mask), as.integer(batch_size))
    } else list(scores = matrix(NA_real_, ncol(kernels), ncol(y)),
      signal_status = rep("insufficient_observations", ncol(y)),
      task_rank = rep(NA_integer_, ncol(kernels)), nuisance_rank = NA_integer_,
      n_observations = sum(mask), comparable = FALSE)
    per_run[[r]] <- result$scores
    design_manifest[[r]] <- data.frame(run = names(bold)[r], candidate = candidate_ids,
      task_rank = result$task_rank, nuisance_rank = result$nuisance_rank,
      n_observations = result$n_observations, comparable = result$comparable)
    event_manifest[[r]] <- data.frame(run = names(bold)[r], voxel = voxel_ids,
      status = result$signal_status)
  }
  # Relative squared shape distance is bounded by two; candidate indices have
  # no metric meaning. Penalizing an index difference would be arbitrary.
  distance <- outer(seq_len(ncol(kernels)), seq_len(ncol(kernels)), Vectorize(function(i, j) {
    sum((kernels[, i] - kernels[, j])^2) / sum(kernels[, i]^2 + kernels[, j]^2)
  }))
  selections <- anchors <- candidate_scores <- vector("list", length(bold))
  training_runs <- vector("list", length(bold))
  names(selections) <- names(anchors) <- names(candidate_scores) <- names(training_runs) <- names(bold)
  for (target in seq_along(bold)) {
    train <- setdiff(seq_along(bold), target)
    training_runs[[target]] <- names(bold)[train]
    # Sum only training matrices, not total-minus-target: even extreme target
    # values cannot enter a fold through cancellation or missingness decisions.
    total <- matrix(0, ncol(kernels), length(voxel_ids))
    counts <- integer(length(voxel_ids))
    for (r in train) {
      usable <- is.finite(per_run[[r]][1L, ])
      total[, usable] <- total[, usable, drop = FALSE] + per_run[[r]][, usable, drop = FALSE]
      counts <- counts + usable
    }
    scores <- sweep(total, 2L, pmax(counts, 1L), "/")
    scores[, counts < min_training_runs] <- NA_real_
    dimnames(scores) <- list(candidate_ids, voxel_ids)
    anchor <- rep(NA_integer_, length(voxel_ids))
    parcel_n <- integer(length(voxel_ids))
    parcel_rows <- list()
    if (!is.null(spatial_groups)) for (g in unique(spatial_groups)) {
      voxels <- which(spatial_groups == g)
      eligible <- voxels[counts[voxels] >= min_training_runs]
      pooled <- if (length(eligible)) rowMeans(scores[, eligible, drop = FALSE]) else rep(NA_real_, ncol(kernels))
      best <- if (all(is.finite(pooled))) which.max(pooled) else NA_integer_
      reliable <- length(eligible) >= min_parcel_voxels && !is.na(best) && pooled[best] >= min_r2
      if (reliable) anchor[voxels] <- best
      parcel_n[voxels] <- length(eligible)
      parcel_rows[[length(parcel_rows) + 1L]] <- data.frame(target_run = names(bold)[target],
        parcel = g, n_voxels = length(eligible), kernel_index = if (reliable) best else fallback_kernel,
        adjusted_partial_r2 = if (!is.na(best)) pooled[best] else NA_real_,
        status = if (reliable) "supported" else "fallback_prespecified")
    }
    anchors[[target]] <- if (length(parcel_rows)) do.call(rbind, parcel_rows) else data.frame()
    chosen <- rep(as.integer(fallback_kernel), length(voxel_ids))
    reason <- rep("fallback_prespecified", length(voxel_ids))
    best <- rep(NA_integer_, length(voxel_ids)); margin <- rep(NA_real_, length(voxel_ids))
    selected_score <- best_score <- rep(NA_real_, length(voxel_ids))
    supported <- which(!is.na(anchor))
    chosen[supported] <- anchor[supported]
    reason[supported] <- if (pooling == "parcel") "pooled_parcel" else "fallback_parcel"
    eligible <- which(counts >= min_training_runs)
    if (length(eligible)) {
      best[eligible] <- max.col(t(scores[, eligible, drop = FALSE]), ties.method = "first")
      best_score[eligible] <- scores[cbind(best[eligible], eligible)]
      if (ncol(kernels) > 1L) {
        second <- rep(-Inf, length(eligible))
        for (k in seq_len(ncol(kernels)))
          second <- pmax(second, ifelse(best[eligible] == k, -Inf, scores[k, eligible]))
        margin[eligible] <- best_score[eligible] - second
      }
      if (pooling == "voxel") {
        objective <- scores[, eligible, drop = FALSE]
        penalized <- which(!is.na(anchor[eligible]))
        if (spatial_strength > 0 && length(penalized))
          objective[, penalized] <- objective[, penalized, drop = FALSE] -
            spatial_strength * distance[, anchor[eligible[penalized]], drop = FALSE]
        candidates <- max.col(t(objective), ties.method = "first")
        keep <- scores[cbind(candidates, eligible)] >= min_r2
        chosen[eligible[keep]] <- candidates[keep]
        reason[eligible[keep]] <- ifelse(spatial_strength > 0 & !is.na(anchor[eligible[keep]]),
          "voxel_spatial", "voxel")
      }
      selected_score[eligible] <- scores[cbind(chosen[eligible], eligible)]
    }
    selections[[target]] <- data.frame(target_run = names(bold)[target], voxel = voxel_ids,
      parcel = if (is.null(spatial_groups)) NA_character_ else spatial_groups,
      kernel_index = chosen, kernel = candidate_ids[chosen], status = reason,
      n_training_runs = counts, n_parcel_voxels = parcel_n, parcel_kernel_index = anchor,
      best_kernel_index = best, best_adjusted_partial_r2 = best_score,
      selected_adjusted_partial_r2 = selected_score, score_margin = margin,
      voxel_evidence = ifelse(counts < min_training_runs, "insufficient_training",
        ifelse(best_score < min_r2, "low_task_support", "supported")))
    if (return_scores) candidate_scores[[target]] <- scores
  }
  structure(list(kernels = kernels, parameters = parameters, selections = selections,
    parcels = anchors, training_runs = training_runs,
    signal_status = do.call(rbind, event_manifest), designs = do.call(rbind, design_manifest),
    scores = if (return_scores) candidate_scores else NULL, run_ids = names(bold),
    voxel_ids = voxel_ids, n_time = vapply(bold, nrow, integer(1)),
    specification = list(events = events, nuisance = nuisance, valid_rows = valid_rows,
      spatial_groups = spatial_groups),
    settings = list(tr = tr, time_offset = offsets, pooling = pooling, amplitude_model = amplitude_model,
      spatial_strength = spatial_strength, fallback_kernel = as.integer(fallback_kernel),
      min_r2 = min_r2, min_training_runs = min_training_runs, min_parcel_voxels = min_parcel_voxels,
      batch_size = as.integer(batch_size), score = "mean_adjusted_partial_r2",
      rank_tolerance = 1e-8,
      nuisance = "run_specific_intercept_linear_trend_and_supplied", noise_model = "OLS",
      selection = "leave_target_run_out", parcel_fallback = !is.null(spatial_groups),
      spatial_regularization = if (pooling == "parcel") "shared_parcel_hrf" else
        if (spatial_strength > 0) "discrete_parcel_shape_penalty" else "none")),
    class = "crossrun_hrf")
}

#' Deconvolve runs using frozen HRFs learned from other runs
#'
#' @param bold Named time-by-voxel run matrices with IDs/order/dimensions matching
#'   the HRF fit. These should be the same preprocessed BOLD used for estimation.
#' @param hrf_fit Object returned by estimate_crossrun_hrf().
#' @param lambda Required fixed nonnegative regularization strength. No target-run
#'   HRF or lambda search is performed. Choose externally or within training data.
#' @param scaling Input scaling: zscore (default), demean, none, or psc.
#' @param penalty Quadratic penalty: diff2, diff1, or ridge.
#' @param backend Fixed linear inverse backend: cpp, r, or scalar wiener.
#' @param batch_size Maximum voxel columns per scaling/multiplication block.
#' @param tail_padding Number of appended zero BOLD observations, default zero.
#'   Use nrow(hrf_fit$kernels) to reproduce the current voxelwise wrapper assumption.
#' @param ridge_floor Numerical diagonal ridge for cpp/r.
#' @param operator_cache_mb Maximum inverse-operator cache size in MiB, default
#'   128. Zero disables reuse across runs. A single operator may exceed the cache
#'   budget while it is being applied, but is not retained in the cache.
#' @return A list with activity matrices, per-run scaling/status tables, the HRF
#'   fit and frozen inverse settings. Activity matrices retain voxel IDs and carry
#'   tr and time_offset attributes for subsequent event alignment.
#' @details One inverse operator is constructed per selected HRF/run length,
#'   cached for consecutive runs of the same length, and applied using block
#'   matrix products for cpp/r. Scalar Wiener uses direct block FFT fits.
#'   HRFs are frozen; no voxelwise factorization or R fitting loop is required.
#'   Nonfinite/constant target signals remain NA. Training observation masks do
#'   not refit a masked inverse; a nonfinite target observation excludes that whole
#'   voxel run. Training nuisance regressors improve HRF selection but are not
#'   automatically regressed out of target BOLD. Preprocess inputs consistently.
#'   Scaling is applied once before optional padding, and its centers/divisors are
#'   returned. Different run z-scores do not preserve absolute between-run units.
#' @export
deconvolve_crossrun_hrf <- function(bold, hrf_fit, lambda, scaling = "zscore",
    penalty = c("diff2", "diff1", "ridge"), backend = c("cpp", "r", "wiener"),
    batch_size = 2048L, tail_padding = 0L, ridge_floor = 1e-8, operator_cache_mb = 128) {
  hrf_validate_bold(bold)
  if (!inherits(hrf_fit, "crossrun_hrf")) stop("hrf_fit must be a crossrun_hrf object")
  if (!identical(names(bold), hrf_fit$run_ids) || !identical(colnames(bold[[1L]]), hrf_fit$voxel_ids) ||
      !identical(vapply(bold, nrow, integer(1)), hrf_fit$n_time))
    stop("BOLD run/voxel IDs and dimensions must match the HRF fit")
  checkmate::assert_number(lambda, lower = 0, finite = TRUE)
  scaling <- match.arg(scaling, c("zscore", "demean", "none", "psc"))
  penalty <- match.arg(penalty); backend <- match.arg(backend)
  checkmate::assert_integerish(batch_size, lower = 1, len = 1L)
  checkmate::assert_integerish(tail_padding, lower = 0, len = 1L)
  checkmate::assert_number(ridge_floor, lower = 0, finite = TRUE)
  checkmate::assert_number(operator_cache_mb, lower = 0, finite = TRUE)
  activity <- scaling_tables <- vector("list", length(bold))
  names(activity) <- names(scaling_tables) <- names(bold)
  operators <- new.env(parent = emptyenv())
  cached_length <- NA_integer_
  cache_keys <- character()
  cache_bytes <- 0
  for (r in seq_along(bold)) {
    y <- bold[[r]]; n <- nrow(y)
    if (!identical(n, cached_length)) {
      rm(list = ls(operators), envir = operators)
      cache_keys <- character(); cache_bytes <- 0
      cached_length <- n
    }
    selected <- hrf_fit$selections[[r]]$kernel_index
    status <- decon_signal_status(y)
    out <- matrix(NA_real_, n, ncol(y), dimnames = dimnames(y))
    centers <- divisors <- rep(NA_real_, ncol(y))
    for (k in unique(selected[status == "ok"])) {
      cache_key <- paste(n, k, sep = ":")
      if (backend != "wiener" && !exists(cache_key, envir = operators, inherits = FALSE)) {
        impulse_basis <- rbind(diag(n), matrix(0, tail_padding, n))
        operator <- deconvolve_reglin(impulse_basis, hrf_fit$kernels[, k], lambda = lambda,
          scaling = "none", backend = backend, penalty = penalty, ridge_floor = ridge_floor)
        operator <- operator[seq_len(n), , drop = FALSE]
        bytes <- as.numeric(object.size(operator))
        if (bytes <= operator_cache_mb * 1024^2) {
          while (length(cache_keys) && cache_bytes + bytes > operator_cache_mb * 1024^2) {
            cache_bytes <- cache_bytes - as.numeric(object.size(get(cache_keys[1], envir = operators)))
            rm(list = cache_keys[1], envir = operators)
            cache_keys <- cache_keys[-1]
          }
          assign(cache_key, operator, envir = operators)
          cache_keys <- c(cache_keys, cache_key); cache_bytes <- cache_bytes + bytes
        }
      } else if (backend != "wiener") {
        operator <- get(cache_key, envir = operators, inherits = FALSE)
      }
      columns <- which(selected == k & status == "ok")
      blocks <- split(columns, ceiling(seq_along(columns) / batch_size))
      for (indices in blocks) {
        scaled <- scale_decon_input(y[, indices, drop = FALSE], scaling = scaling)
        out[, indices] <- if (backend == "wiener") {
          # Preserve the FFT backend's O(T log T) application cost instead of
          # replacing it with a dense O(T^2) operator multiplication.
          fit_input <- rbind(scaled$y, matrix(0, tail_padding, length(indices)))
          deconvolve_reglin(fit_input, hrf_fit$kernels[, k], lambda = lambda,
            scaling = "none", backend = "wiener")[seq_len(n), , drop = FALSE]
        } else operator %*% scaled$y
        centers[indices] <- scaled$center; divisors[indices] <- scaled$scale
      }
    }
    attr(out, "tr") <- hrf_fit$settings$tr
    attr(out, "time_offset") <- hrf_fit$settings$time_offset[r]
    activity[[r]] <- out
    scaling_tables[[r]] <- data.frame(run = names(bold)[r], voxel = colnames(y),
      status = status, center = centers, scale = divisors, scaling = scaling)
  }
  list(activity = activity, scaling = scaling_tables, hrf_fit = hrf_fit,
    settings = list(lambda = lambda, scaling = scaling, penalty = penalty,
      backend = backend, tail_padding = as.integer(tail_padding), ridge_floor = ridge_floor,
      batch_size = as.integer(batch_size), operator_cache_mb = operator_cache_mb,
      hrf_selection = "frozen_crossrun"))
}
