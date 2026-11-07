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
#' @param noise_model Default ar1 uses pooled, grid-REML covariance estimates;
#'   ols disables whitening. Pooling uses supplied spatial_groups or the run.
#' @param noise_rho Optional fixed AR(1) coefficient, scalar or named in run
#'   order. Bypasses noise estimation. Requires noise_model = "ar1".
#' @param noise_grid Ordered AR(1) grid within [-0.95, 0.95], default -0.8 to 0.9
#'   in steps of 0.1. Boundary estimates are explicitly labeled.
#' @param evidence Default calibrated uses conditional GLS F tests with library
#'   corrections, prespecified-shape comparison, and training-run compatibility.
#'   Legacy retains the adjusted-R-squared heuristic. Conditional tests are exact
#'   for known Gaussian AR covariance and a correct mean model; estimated noise,
#'   non-Gaussian errors, and model misspecification require empirical validation.
#' @param evidence_alpha Per-voxel evidence threshold, default 0.05. Corrects the
#'   candidate search and compatibility comparisons across training runs, not
#'   multiple testing across the brain. Not a posterior probability.
#' @param canonical_strength Fixed shape-distance penalty toward fallback_kernel
#'   in calibrated mode, default 0.01. Zero removes this additional shrinkage.
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
#'   Calibrated mode additionally penalizes distance from the prespecified shape.
#'   The legacy gate accepted null free-event fits in colored-noise simulations.
#'   Calibrated mode tests task support, departures from the prespecified shape
#'   using pairwise encompassing models, and compatibility with every usable
#'   training run. Fisher combination uses independent runs; Bonferroni corrections
#'   account for candidate comparisons. Unsupported or inconsistent shapes fall
#'   back. Parcel anchors require min_parcel_voxels individually supported voxels;
#'   voxel p-values are not combined as independent spatial observations.
#'   Run compatibility cannot be assessed with a single training run. There is
#'   no access to target-run vascular shape or guarantee of transfer to that run.
#'   AR covariance is estimated per run using a fixed event-amplitude mean model
#'   and, when residual degrees of freedom allow, its first two temporal
#'   derivatives. Noise likelihood profiles are pooled, not the BOLD signals.
#'   All whitening uses the original acquisition clock, including gaps. No
#'   automatic nested tuning or physiological amplitude calibration is included.
#'   Nonfinite/constant/nuisance-only voxel runs are omitted; observation masks are
#'   common across voxels. Candidate task ranks must match within a run and leave
#'   residual degrees of freedom, otherwise that run cannot inform HRF selection.
#' @export
estimate_crossrun_hrf <- function(bold, events, tr, kernel_grid, nuisance = NULL,
    valid_rows = NULL, time_offset = 0, spatial_groups = NULL,
    pooling = c("voxel", "parcel"), spatial_strength = 0, fallback_kernel = 1L,
    min_r2 = 0.01, min_training_runs = 1L, min_parcel_voxels = 5L,
    batch_size = 2048L, return_scores = FALSE, amplitude_model = c("condition", "trial"),
    noise_model = c("ar1", "ols"), noise_rho = NULL,
    noise_grid = seq(-0.8, 0.9, by = 0.1), evidence = c("calibrated", "legacy"),
    evidence_alpha = 0.05, canonical_strength = 0.01) {
  hrf_validate_bold(bold)
  pooling <- match.arg(pooling)
  amplitude_model <- match.arg(amplitude_model)
  noise_model <- match.arg(noise_model)
  evidence <- match.arg(evidence)
  checkmate::assert_number(evidence_alpha, lower = .Machine$double.eps, upper = 0.5, finite = TRUE)
  checkmate::assert_number(canonical_strength, lower = 0, finite = TRUE)
  checkmate::assert_numeric(noise_grid, lower = -0.95, upper = 0.95, min.len = 2L,
    any.missing = FALSE, finite = TRUE, unique = TRUE, sorted = TRUE)
  if (noise_model == "ols" && !is.null(noise_rho)) stop("noise_rho requires noise_model = 'ar1'")
  if (!is.null(noise_rho)) {
    noise_rho <- hrf_run_numeric(noise_rho, bold, "noise_rho")
    checkmate::assert_numeric(noise_rho, lower = -0.95, upper = 0.95)
  }
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
  per_run <- run_results <- noise <- design_manifest <- event_manifest <- vector("list", length(bold))
  names(per_run) <- names(run_results) <- names(noise) <- names(design_manifest) <- names(event_manifest) <- names(bold)
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
    if (any(!is.finite(z[mask, , drop = FALSE]))) stop("Nuisance values must be finite on usable rows")
    noise[[r]] <- hrf_estimate_noise(y, ev, tr, offsets[r], kernels, fallback_kernel,
      z, mask, spatial_groups, noise_grid, batch_size,
      if (noise_model == "ols") 0 else if (!is.null(noise_rho)) noise_rho[r] else NULL)
    result <- list(scores = matrix(NA_real_, ncol(kernels), ncol(y)),
      signal_status = rep("insufficient_observations", ncol(y)))
    if (evidence == "calibrated") for (field in c("task_p", "fallback_p", "incompatible_p"))
      result[[field]] <- matrix(NA_real_, ncol(kernels), ncol(y))
    manifest <- list()
    for (rho in unique(noise[[r]]$rho)) {
      indices <- which(noise[[r]]$rho == rho)
      part <- if (sum(mask) >= 3L) {
        hrf_score_run_cpp(y[, indices, drop = FALSE], design, z, which(mask),
          as.integer(batch_size), rho = rho, evidence = evidence == "calibrated",
          fallback_index = as.integer(fallback_kernel))
      } else list(scores = matrix(NA_real_, ncol(kernels), length(indices)),
        signal_status = rep("insufficient_observations", length(indices)),
        task_rank = rep(NA_integer_, ncol(kernels)), nuisance_rank = NA_integer_,
        n_observations = sum(mask), comparable = FALSE)
      result$scores[, indices] <- part$scores
      result$signal_status[indices] <- part$signal_status
      if (evidence == "calibrated" && sum(mask) >= 3L)
        for (field in c("task_p", "fallback_p", "incompatible_p")) result[[field]][, indices] <- part[[field]]
      manifest[[length(manifest) + 1L]] <- data.frame(run = names(bold)[r], candidate = candidate_ids,
        noise_rho = rho, task_rank = part$task_rank, nuisance_rank = part$nuisance_rank,
        n_observations = part$n_observations, comparable = part$comparable)
    }
    if (noise_model == "ar1" && evidence == "calibrated" && sum(mask) >= 3L) {
      unknown_noise <- grepl("^insufficient_noise", noise[[r]]$method) &
        !result$signal_status %in% c("nonfinite", "constant", "nuisance_only")
      result$signal_status[unknown_noise] <- "noise_unidentified"
      for (field in c("scores", "task_p", "fallback_p", "incompatible_p"))
        result[[field]][, unknown_noise] <- NA_real_
    }
    run_results[[r]] <- result
    per_run[[r]] <- result$scores
    design_manifest[[r]] <- do.call(rbind, manifest)
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
    gate <- if (evidence == "calibrated") hrf_fold_evidence(run_results, train, counts,
      scores, fallback_kernel, evidence_alpha, min_r2) else NULL
    anchor <- rep(NA_integer_, length(voxel_ids))
    parcel_n <- integer(length(voxel_ids))
    parcel_rows <- list()
    if (!is.null(spatial_groups)) for (g in unique(spatial_groups)) {
      voxels <- which(spatial_groups == g)
      eligible <- voxels[counts[voxels] >= min_training_runs]
      pooled <- if (length(eligible)) rowMeans(scores[, eligible, drop = FALSE]) else rep(NA_real_, ncol(kernels))
      admissible <- if (is.null(gate)) rep(TRUE, ncol(kernels)) else
        rowSums(gate$supported[, eligible, drop = FALSE]) >= min_parcel_voxels
      pooled_objective <- pooled
      if (!is.null(gate)) pooled_objective <- pooled_objective - canonical_strength * distance[, fallback_kernel]
      pooled_objective[!admissible] <- -Inf
      best <- if (any(is.finite(pooled_objective))) which.max(pooled_objective) else NA_integer_
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
    if (!is.null(gate) && length(supported)) supported <- supported[
      gate$compatible[cbind(anchor[supported], supported)]]
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
        if (!is.null(gate)) {
          objective <- objective - canonical_strength * distance[, fallback_kernel]
          objective[!gate$supported[, eligible, drop = FALSE]] <- -Inf
        }
        penalized <- which(!is.na(anchor[eligible]))
        if (spatial_strength > 0 && length(penalized))
          objective[, penalized] <- objective[, penalized, drop = FALSE] -
            spatial_strength * distance[, anchor[eligible[penalized]], drop = FALSE]
        candidates <- max.col(t(objective), ties.method = "first")
        keep <- scores[cbind(candidates, eligible)] >= min_r2
        if (!is.null(gate)) keep <- keep & gate$supported[cbind(candidates, eligible)]
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
    if (!is.null(gate)) {
      sel <- selections[[target]]
      idx <- cbind(sel$kernel_index, seq_along(voxel_ids))
      sel$selected_task_p <- gate$task_p[idx]
      sel$selected_shape_p <- gate$shape_p[idx]
      sel$n_compatible_shapes <- colSums(gate$compatible)
      sel$n_supported_shapes <- colSums(gate$supported)
      sel$stability <- ifelse(counts < 2L, "not_assessable", ifelse(
        sel$n_compatible_shapes == 0L, "no_common_library_shape", "compatible_training_runs"))
      own_support <- gate$supported[idx]
      sel$voxel_evidence <- ifelse(counts < min_training_runs, "insufficient_training",
        ifelse(sel$n_compatible_shapes == 0L, "inconsistent_or_misspecified",
          ifelse(!own_support, "unidentified",
            ifelse(sel$kernel_index == fallback_kernel, "prespecified_shape", "supported_departure"))))
      selections[[target]] <- sel
    }
    if (return_scores) candidate_scores[[target]] <- scores
  }
  structure(list(kernels = kernels, parameters = parameters, selections = selections,
    parcels = anchors, training_runs = training_runs,
    signal_status = do.call(rbind, event_manifest), designs = do.call(rbind, design_manifest),
    scores = if (return_scores) candidate_scores else NULL,
    noise = noise, run_evidence = if (return_scores && evidence == "calibrated") run_results else NULL,
    run_ids = names(bold),
    voxel_ids = voxel_ids, n_time = vapply(bold, nrow, integer(1)),
    specification = list(events = events, nuisance = nuisance, valid_rows = valid_rows,
      spatial_groups = spatial_groups),
    settings = list(tr = tr, time_offset = offsets, pooling = pooling, amplitude_model = amplitude_model,
      spatial_strength = spatial_strength, fallback_kernel = as.integer(fallback_kernel),
      min_r2 = min_r2, min_training_runs = min_training_runs, min_parcel_voxels = min_parcel_voxels,
      batch_size = as.integer(batch_size), score = "mean_adjusted_partial_r2",
      rank_tolerance = 1e-8, evidence = evidence, evidence_alpha = evidence_alpha,
      canonical_strength = canonical_strength, noise_grid = noise_grid,
      noise_pooling = if (is.null(spatial_groups)) "run" else "supplied_parcel",
      nuisance = "run_specific_intercept_linear_trend_and_supplied", noise_model = noise_model,
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
#' @param nuisance_mode Jointly fit unpenalized intercept, trend, and supplied
#'   nuisance terms, or ignore them (default, preserving prepared input signals).
#'   Joint projection is applied to both data and design, but can amplify weakly
#'   identified low-frequency activity in a flexible continuous inverse. Validate
#'   it for the intended design rather than assuming it improves recovery.
#' @param observation_mask Use fit's valid_rows in the inverse (fit), or all
#'   observations (all). Defaults to fit for AR/calibrated, all for OLS/legacy.
#' @return A list with activity matrices, per-run scaling/status tables, the HRF
#'   fit and frozen inverse settings. Activity matrices retain voxel IDs and carry
#'   tr and time_offset attributes for subsequent event alignment.
#' @details One inverse operator is constructed per selected HRF/run length,
#'   cached for consecutive runs with identical masks/nuisance designs, and applied using block
#'   matrix products for cpp/r. Scalar Wiener uses direct block FFT fits.
#'   HRFs are frozen; no voxelwise factorization or R fitting loop is required.
#'   AR/calibrated fits use a masked GLS inverse with the frozen run-specific
#'   noise covariance. A gap of d acquisitions uses rho^d and the corresponding
#'   accumulated innovation variance. Censored rows contribute no observations;
#'   latent activity is still estimated on the complete original grid. Its
#'   observation_mask attribute identifies modeled versus observed samples.
#'   Nonfinite/constant signals on usable observations remain NA. Nonfinite
#'   censored observations do not invalidate a run. Joint nuisance fitting avoids
#'   projecting only the BOLD and inadvertently changing the task model.
#'   Wiener requires OLS, nuisance_mode = "ignore", observation_mask = "all".
#'   Scaling is applied once before optional padding, and its centers/divisors are
#'   returned. Different run z-scores do not preserve absolute between-run units.
#' @export
deconvolve_crossrun_hrf <- function(bold, hrf_fit, lambda, scaling = "zscore",
    penalty = c("diff2", "diff1", "ridge"), backend = c("cpp", "r", "wiener"),
    batch_size = 2048L, tail_padding = 0L, ridge_floor = 1e-8, operator_cache_mb = 128,
    nuisance_mode = NULL, observation_mask = NULL) {
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
  modern <- identical(hrf_fit$settings$noise_model, "ar1") || identical(hrf_fit$settings$evidence, "calibrated")
  if (is.null(nuisance_mode)) nuisance_mode <- "ignore"
  if (is.null(observation_mask)) observation_mask <- if (modern) "fit" else "all"
  nuisance_mode <- match.arg(nuisance_mode, c("joint", "ignore"))
  observation_mask <- match.arg(observation_mask, c("fit", "all"))
  if (backend == "wiener" && (identical(hrf_fit$settings$noise_model, "ar1") ||
      nuisance_mode != "ignore" || observation_mask != "all"))
    stop("Wiener does not implement masked, nuisance-adjusted AR GLS; use cpp/r, or an OLS fit with nuisance_mode = 'ignore', observation_mask = 'all'")
  settings <- list(lambda = lambda, scaling = scaling, penalty = penalty,
    backend = backend, tail_padding = as.integer(tail_padding), ridge_floor = ridge_floor,
    batch_size = as.integer(batch_size), operator_cache_mb = operator_cache_mb,
    hrf_selection = "frozen_crossrun", nuisance_mode = nuisance_mode,
    observation_mask = observation_mask, noise_model = hrf_fit$settings$noise_model)
  activity <- scaling_tables <- vector("list", length(bold))
  names(activity) <- names(scaling_tables) <- names(bold)
  operators <- new.env(parent = emptyenv())
  cached_design <- NULL
  cache_keys <- character(); cache_bytes <- 0
  for (r in seq_along(bold)) {
    y <- bold[[r]]; n <- nrow(y)
    mask <- if (observation_mask == "fit") hrf_fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    z <- if (nuisance_mode == "joint") cbind(1, seq(-1, 1, length.out = n),
      hrf_fit$specification$nuisance[[r]]) else matrix(numeric(), n, 0L)
    design <- list(n = n, rows = rows, nuisance = z)
    if (!identical(design, cached_design)) {
      rm(list = ls(operators), envir = operators)
      cache_keys <- character(); cache_bytes <- 0
      cached_design <- design
    }
    selected <- hrf_fit$selections[[r]]$kernel_index
    rho <- if (is.null(hrf_fit$noise)) rep(0, ncol(y)) else hrf_fit$noise[[r]]$rho
    status <- if (length(rows) >= 3L) decon_signal_status(y[rows, , drop = FALSE]) else
      rep("insufficient_observations", ncol(y))
    if (identical(hrf_fit$settings$noise_model, "ar1") && !is.null(hrf_fit$noise)) {
      unknown <- grepl("^insufficient_noise", hrf_fit$noise[[r]]$method)
      status[status == "ok" & unknown] <- "noise_unidentified"
    }
    out <- matrix(NA_real_, n, ncol(y), dimnames = dimnames(y))
    centers <- divisors <- rep(NA_real_, ncol(y))
    strata <- paste(selected, sprintf("%.17g", rho), sep = ":")
    for (key in unique(strata[status == "ok"])) {
      columns <- which(strata == key & status == "ok")
      k <- selected[columns[1L]]; ar <- rho[columns[1L]]
      if (backend != "wiener") {
        if (exists(key, envir = operators, inherits = FALSE)) operator <- get(key, envir = operators) else {
          operator <- hrf_inverse_operator(n, hrf_fit$kernels[, k], rows, ar, z, settings)
          bytes <- as.numeric(object.size(operator))
          if (bytes <= operator_cache_mb * 1024^2) {
            while (length(cache_keys) && cache_bytes + bytes > operator_cache_mb * 1024^2) {
              cache_bytes <- cache_bytes - as.numeric(object.size(get(cache_keys[1L], envir = operators)))
              rm(list = cache_keys[1L], envir = operators); cache_keys <- cache_keys[-1L]
            }
            assign(key, operator, envir = operators)
            cache_keys <- c(cache_keys, key); cache_bytes <- cache_bytes + bytes
          }
        }
      }
      for (indices in split(columns, ceiling(seq_along(columns) / batch_size))) {
        scaled <- scale_decon_input(y[rows, indices, drop = FALSE], scaling = scaling)
        out[, indices] <- if (backend == "wiener") {
          fit_input <- rbind(scaled$y, matrix(0, tail_padding, length(indices)))
          deconvolve_reglin(fit_input, hrf_fit$kernels[, k], lambda = lambda,
            scaling = "none", backend = "wiener")[seq_len(n), , drop = FALSE]
        } else operator %*% scaled$y
        centers[indices] <- scaled$center; divisors[indices] <- scaled$scale
      }
    }
    attr(out, "tr") <- hrf_fit$settings$tr
    attr(out, "time_offset") <- hrf_fit$settings$time_offset[r]
    attr(out, "observation_mask") <- seq_len(n) %in% rows
    activity[[r]] <- out
    scaling_tables[[r]] <- data.frame(run = names(bold)[r], voxel = colnames(y),
      status = status, center = centers, scale = divisors, scaling = scaling,
      noise_rho = rho, n_observations = length(rows))
  }
  list(activity = activity, scaling = scaling_tables, hrf_fit = hrf_fit, settings = settings)
}
