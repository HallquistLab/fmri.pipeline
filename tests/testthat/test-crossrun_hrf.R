# Freeze the original OLS/heuristic contract; new defaults have separate tests.
legacy_estimate_crossrun_hrf <- function(...) {
  estimate_crossrun_hrf(..., noise_model = "ols", evidence = "legacy")
}

# Independent, small-grid forward reference; no native design/scoring helpers.
hrf_reference_design <- function(n, tr, offset, events, kernels) {
  history <- nrow(kernels) - 1L
  neural <- matrix(0, n + history, length(unique(events$condition)))
  for (e in seq_len(nrow(events))) {
    position <- (events$onset[e] - offset) / tr + history
    duration <- if (is.null(events$duration)) 0 else events$duration[e]
    j <- match(events$condition[e], unique(events$condition))
    for (t in 0:(n + history - 1L)) {
      value <- if (duration == 0) max(0, 1 - abs(t - position)) else {
        max(0, min(position + duration / tr, t + 1) - max(position, t))
      }
      neural[t + 1L, j] <- neural[t + 1L, j] + value
    }
  }
  lapply(seq_len(ncol(kernels)), function(k) {
    vapply(seq_len(ncol(neural)), function(j) {
      convolved <- convolve(neural[, j], rev(kernels[, k]), type = "open")
      convolved[history + seq_len(n)]
    }, numeric(n))
  })
}

hrf_test_data <- function() {
  set.seed(619)
  kernels <- cbind(early = c(0, .15, .5, .3, .05, 0, 0),
    late = c(0, 0, .05, .15, .5, .3, 0))
  events <- list(A = data.frame(onset = c(-1.5, 3, 12.5, 21, 32, 45, 56.5, 68, 80),
    condition = rep(c("cue", "feedback", "response"), 3)))
  events$B <- transform(events$A, onset = onset + c(0, 1, 2, 0, 1, 0, 2, 1, 0))
  events$C <- transform(events$B, onset = onset + 0.5)
  truth <- rep(1:2, each = 4)
  bold <- lapply(events, function(ev) {
    designs <- hrf_reference_design(96, 1, 0, ev, kernels)
    y <- vapply(seq_along(truth), function(v) {
      as.numeric(designs[[truth[v]]] %*% c(1, 2, -0.5)) * (1 + v / 4) +
        100 + seq(-1, 1, length.out = 96) * v + rnorm(96, sd = .002)
    }, numeric(96))
    colnames(y) <- paste0("v", 1:8)
    y
  })
  groups <- stats::setNames(rep(c("p1", "p2"), each = 4), colnames(bold$A))
  list(bold = bold, events = events, kernels = kernels, truth = truth, groups = groups)
}

test_that("native event design agrees with causal convolution including pre-run history", {
  kernels <- cbind(c(0, .4, .6, 0), c(0, .2, .3, .5))
  ev <- data.frame(onset = c(-3, -0.5, 3.75, 9, 100),
    condition = c("a", "b", "a", "b", "a"), duration = c(0, 2.5, 0, 3.2, 0))
  expected <- hrf_reference_design(20, .75, 1.5, ev, kernels)
  got <- fmri.pipeline:::hrf_event_design_cpp(20L, .75, 1.5, ev$onset,
    ev$duration, match(ev$condition, unique(ev$condition)), 2L, kernels)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("native scores match independent masked nuisance-adjusted least squares", {
  d <- hrf_test_data()
  y <- d$bold$A
  nuisance <- cbind(1, seq(-1, 1, length.out = 96), sin(1:96 / 9))
  mask <- rep(TRUE, 96); mask[c(2, 20:23, 70)] <- FALSE
  y[!mask, 1] <- NA_real_
  y[10, 2] <- NA_real_; y[, 3] <- 10; y[, 4] <- nuisance[, 2]
  design <- hrf_reference_design(96, 1, 0, d$events$A, d$kernels)
  got <- fmri.pipeline:::hrf_score_run_cpp(y, design, nuisance, which(mask), 3L)
  for (k in 1:2) for (v in c(1, 5:8)) {
    base <- lm.fit(nuisance[mask, ], y[mask, v])
    full <- lm.fit(cbind(nuisance[mask, ], design[[k]][mask, ]), y[mask, v])
    expected <- 1 - (sum(full$residuals^2) / (sum(mask) - full$rank)) /
      (sum(base$residuals^2) / (sum(mask) - base$rank))
    expect_equal(got$scores[k, v], expected, tolerance = 1e-9)
  }
  expect_equal(got$signal_status[2:4], c("nonfinite", "constant", "nuisance_only"))
  expect_true(all(is.na(got$scores[, 2:4])))
  again <- fmri.pipeline:::hrf_score_run_cpp(y, design, nuisance, which(mask), 100L)
  expect_equal(again, got, tolerance = 1e-12)
  scaled_z <- sweep(nuisance, 2, c(1e-10, 1e10, 1e-5), "*")
  rescaled <- fmri.pipeline:::hrf_score_run_cpp(y, design, scaled_z, which(mask), 3L)
  expect_equal(rescaled, got, tolerance = 1e-10)
})

test_that("cross-run shapes and parcel pooling recover known HRFs", {
  d <- hrf_test_data()
  fit <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels,
    spatial_groups = d$groups, pooling = "parcel", min_parcel_voxels = 2,
    return_scores = TRUE, batch_size = 3)
  for (r in names(d$bold)) {
    expect_equal(fit$selections[[r]]$kernel_index, d$truth)
    expect_true(all(fit$selections[[r]]$status == "pooled_parcel"))
    expect_true(all(fit$selections[[r]]$n_training_runs == 2L))
    expect_false(r %in% fit$training_runs[[r]])
    expect_true(all(fit$parcels[[r]]$status == "supported"))
    expect_true(all(fit$selections[[r]]$selected_adjusted_partial_r2 > .99))
  }
  expect_equal(colSums(fit$kernels), c(early = 1, late = 1))
})

test_that("target BOLD, events, mask and nuisance do not affect its HRF selection", {
  d <- hrf_test_data()
  z <- lapply(d$bold, function(y) matrix(sin(1:96), 96, 1))
  masks <- lapply(d$bold, function(y) rep(TRUE, 96))
  settings <- list(bold = d$bold, events = d$events, tr = 1, kernel_grid = d$kernels,
    spatial_groups = d$groups, spatial_strength = .1, min_parcel_voxels = 2,
    nuisance = z, valid_rows = masks, return_scores = TRUE)
  original <- do.call(legacy_estimate_crossrun_hrf, settings)
  settings$bold$A[,] <- NA_real_
  settings$events$A$onset <- settings$events$A$onset + 1000
  settings$valid_rows$A[] <- FALSE
  settings$nuisance$A[,] <- cos(1:96)
  changed <- do.call(legacy_estimate_crossrun_hrf, settings)
  expect_identical(changed$selections$A, original$selections$A)
  expect_identical(changed$parcels$A, original$parcels$A)
  expect_identical(changed$scores$A, original$scores$A)
  expect_true(all(changed$signal_status$status[changed$signal_status$run == "A"] == "insufficient_observations"))
})

test_that("spatial shape penalty follows its declared objective and does not mix parcels", {
  d <- hrf_test_data()
  # One voxel disagrees with its parcel, but has good independent task support.
  for (r in names(d$bold)) d$bold[[r]][, 4] <- d$bold[[r]][, 5]
  plain <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels,
    spatial_groups = d$groups, min_parcel_voxels = 2, return_scores = TRUE)
  smooth <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels,
    spatial_groups = d$groups, spatial_strength = 100, min_parcel_voxels = 2,
    min_r2 = 0, return_scores = TRUE)
  expect_equal(plain$selections$A$kernel_index[4], 2L)
  for (v in 1:8) {
    prior <- smooth$selections$A$parcel_kernel_index[v]
    distance <- vapply(1:2, function(k) sum((d$kernels[, k] - d$kernels[, prior])^2) /
      sum(d$kernels[, k]^2 + d$kernels[, prior]^2), numeric(1))
    expected <- which.max(smooth$scores$A[, v] - 100 * distance)
    expect_equal(smooth$selections$A$kernel_index[v], unname(expected))
  }
  expect_equal(smooth$selections$A$kernel_index, d$truth)
  expect_true(all(smooth$selections$A$status == "voxel_spatial"))
})

test_that("weak and missing training signals use explicit fallback hierarchy", {
  d <- hrf_test_data()
  for (r in names(d$bold)) d$bold[[r]][, 1] <- NA_real_
  fit <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels,
    spatial_groups = d$groups, min_parcel_voxels = 2)
  expect_equal(fit$selections$A$status[1], "fallback_parcel")
  expect_equal(fit$selections$A$n_training_runs[1], 0L)
  expect_equal(fit$selections$A$kernel_index[1], 1L)
  fallback <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels, fallback_kernel = 2)
  expect_equal(fallback$selections$A$status[1], "fallback_prespecified")
  expect_equal(fallback$selections$A$kernel_index[1], 2L)
  # No task support: all events are far outside the measured run.
  outside <- lapply(d$events, function(ev) transform(ev, onset = onset + 1000))
  nofit <- legacy_estimate_crossrun_hrf(d$bold, outside, 1, d$kernels, fallback_kernel = 2)
  expect_true(all(nofit$selections$A$status == "fallback_prespecified"))
  expect_true(all(nofit$signal_status$status == "unidentifiable_design"))
  expect_true(all(!nofit$designs$comparable))
})

test_that("frozen grouped application matches direct inverses including scaling and padding", {
  d <- hrf_test_data()
  d$bold$A[3, 1] <- NA_real_; d$bold$B[, 2] <- 100
  fit <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels)
  for (backend in c("cpp", "r", "wiener")) for (scaling in c("none", "demean", "zscore", "psc")) {
    got <- deconvolve_crossrun_hrf(d$bold, fit, lambda = .03, scaling = scaling,
      backend = backend, tail_padding = 7L, batch_size = 2L)
    for (r in names(d$bold)) for (k in 1:2) {
      indices <- which(fit$selections[[r]]$kernel_index == k & got$scaling[[r]]$status == "ok")
      if (!length(indices)) next
      scaled <- fmri.pipeline:::scale_decon_input(d$bold[[r]][, indices, drop = FALSE], scaling = scaling)
      expected <- deconvolve_reglin(rbind(scaled$y, matrix(0, 7, length(indices))),
        d$kernels[, k], lambda = .03, scaling = "none", backend = backend)[1:96, , drop = FALSE]
      expect_equal(unname(got$activity[[r]][, indices, drop = FALSE]), unname(expected),
        ignore_attr = TRUE, tolerance = 1e-8)
      expect_equal(got$scaling[[r]]$center[indices], scaled$center, ignore_attr = TRUE)
      expect_equal(got$scaling[[r]]$scale[indices], scaled$scale, ignore_attr = TRUE)
    }
    expect_true(all(is.na(got$activity$A[, 1])))
    expect_true(all(is.na(got$activity$B[, 2])))
    expect_equal(colnames(got$activity$A), colnames(d$bold$A))
    expect_equal(attr(got$activity$A, "tr"), 1)
  }
})

test_that("HRF library and identity guards prevent ambiguous application", {
  kernels <- hrf_library(.8)
  expect_equal(ncol(kernels), 9)
  expect_equal(colSums(kernels), rep(1, 9), ignore_attr = TRUE)
  expect_equal(as.numeric(kernels[, 1]), as.numeric(fmri.pipeline:::spm_hrf(.8)$hrf))
  expect_equal(attr(kernels, "parameters")[1, ],
    data.frame(peak_delay = 6, dispersion = 1, undershoot_ratio = 6))
  d <- hrf_test_data()
  mismatch <- d$bold; colnames(mismatch$B) <- rev(colnames(mismatch$B))
  expect_error(legacy_estimate_crossrun_hrf(mismatch, d$events, 1, d$kernels), "identical voxel IDs")
  expect_error(legacy_estimate_crossrun_hrf(d$bold, unname(d$events), 1, d$kernels), "same run names")
  expect_error(legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels, pooling = "parcel"), "explicit spatial_groups")
  expect_error(legacy_estimate_crossrun_hrf(d$bold, d$events, 1, hrf_library(.8)), "sampling TR")
  fit <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels)
  expect_error(deconvolve_crossrun_hrf(d$bold, fit), "lambda")
  expect_error(deconvolve_crossrun_hrf(mismatch, fit, .1), "identical voxel IDs")
  expect_error(deconvolve_crossrun_hrf(rev(d$bold), fit, .1), "must match")
})

test_that("two runs, unequal lengths and bounded operator caching are supported", {
  d <- hrf_test_data()
  d$bold <- d$bold[1:2]; d$events <- d$events[1:2]
  d$bold$B <- d$bold$B[1:90, , drop = FALSE]
  fit <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels)
  expect_equal(fit$training_runs, list(A = "B", B = "A"))
  expect_true(all(fit$selections$A$n_training_runs == 1))
  no_cache <- deconvolve_crossrun_hrf(d$bold, fit, .1, operator_cache_mb = 0)
  evict <- deconvolve_crossrun_hrf(d$bold, fit, .1, operator_cache_mb = .08)
  expect_equal(no_cache$activity, evict$activity)
  expect_equal(dim(no_cache$activity$B), c(90L, 8L))
  for (r in names(d$bold)) {
    reference <- vapply(1:8, function(v) {
      k <- fit$selections[[r]]$kernel_index[v]
      deconvolve_reglin(d$bold[[r]][, v], d$kernels[, k], lambda = .1)[, 1]
    }, numeric(nrow(d$bold[[r]])))
    expect_equal(as.vector(no_cache$activity[[r]]), as.vector(reference), tolerance = 1e-8)
  }
})

test_that("masks preserve acquisition timing and require comparable candidate ranks", {
  d <- hrf_test_data()
  masks <- lapply(d$bold, function(y) rep(TRUE, nrow(y)))
  masks$B[10:15] <- FALSE; masks$B[50] <- NA
  d$bold$B[10:15, ] <- NA_real_
  fit <- legacy_estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels, valid_rows = masks)
  expect_true(all(fit$signal_status$status == "ok"))
  expect_equal(unique(fit$designs$n_observations[fit$designs$run == "B"]), 89L)
  expect_equal(fit$specification$valid_rows, masks)
  # The second condition is exactly absorbed by nuisance; it must not acquire
  # spurious rank when the residual is normalized after projection.
  x <- cbind(a = sin(1:30), b = cos(1:30))
  y <- matrix(sin(1:30) + rnorm(30, sd = .1), 30, 1)
  result <- fmri.pipeline:::hrf_score_run_cpp(y, list(x, x), cbind(1, x[, 2]), 1:30)
  expect_equal(result$task_rank, c(1L, 1L))
  uneven <- fmri.pipeline:::hrf_score_run_cpp(y, list(x, cbind(x[, 2], x[, 2])),
    cbind(1, x[, 2]), 1:30)
  expect_false(uneven$comparable)
  expect_true(all(is.na(uneven$scores)))
})

test_that("free trial amplitudes identify an HRF when signed responses average away", {
  d <- hrf_test_data()
  ev <- data.frame(onset = seq(10, 145, by = 15), condition = "feedback")
  trial_events <- transform(ev, condition = seq_len(nrow(ev)))
  x <- hrf_reference_design(160, 1, 0, trial_events, d$kernels)
  y <- vapply(1:2, function(k) as.numeric(x[[k]] %*% rep(c(-1, 1), 5)) + 100, numeric(160))
  colnames(y) <- c("early", "late")
  runs <- list(A = y, B = y)
  events <- list(A = ev, B = ev)
  condition <- legacy_estimate_crossrun_hrf(runs, events, 1, d$kernels)
  trial <- legacy_estimate_crossrun_hrf(runs, events, 1, d$kernels,
    amplitude_model = "trial", return_scores = TRUE)
  expect_true(all(condition$selections$A$status == "fallback_prespecified"))
  expect_equal(trial$selections$A$kernel_index, 1:2)
  expect_true(all(trial$selections$A$best_adjusted_partial_r2 > .9999))
  expect_true(all(trial$designs$task_rank == nrow(ev)))
  expect_equal(trial$settings$amplitude_model, "trial")
})

test_that("hrf_library presets cover extended and dense physiological grids", {
  std <- hrf_library(1, preset = "standard")
  ext <- hrf_library(1, preset = "extended")
  den <- hrf_library(1, preset = "dense")
  expect_equal(ncol(std), 9L)
  expect_equal(ncol(ext), 15L)
  expect_equal(ncol(den), 40L)
  expect_equal(attr(ext, "parameters")[1L, ],
    data.frame(peak_delay = 6, dispersion = 1, undershoot_ratio = 6))
  expect_equal(attr(den, "parameters")[1L, ],
    data.frame(peak_delay = 6, dispersion = 1, undershoot_ratio = 6))
  expect_equal(colSums(ext), rep(1, ncol(ext)), ignore_attr = TRUE)
  expect_equal(colSums(den), rep(1, ncol(den)), ignore_attr = TRUE)
})
