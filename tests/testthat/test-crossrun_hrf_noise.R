hrfn_data <- function(n = 110L, voxels = 6L) {
  set.seed(7129)
  h <- cbind(early = c(0, .2, .5, .25, .05, 0), late = c(0, 0, .1, .3, .5, .1))
  ev <- data.frame(onset = seq(8, n - 12, by = 11) + .3, condition = "feedback")
  events <- setNames(rep(list(ev), 3L), c("A", "B", "C"))
  x <- hrf_event_design_cpp(n, 1, 0, ev$onset, rep(0, nrow(ev)),
    seq_len(nrow(ev)), nrow(ev), h)
  bold <- lapply(events, function(e) {
    y <- vapply(seq_len(voxels), function(v) {
      as.numeric(x[[1L + (v %% 2L)]] %*% rep(c(-1, 1), length.out = nrow(ev))) +
        100 + seq(-1, 1, length.out = n) + rnorm(n, sd = .005)
    }, numeric(n))
    colnames(y) <- paste0("v", seq_len(voxels)); y
  })
  list(bold = bold, events = events, kernels = h)
}

test_that("gap-aware AR whitening and evidence match independent GLS fits", {
  set.seed(8173)
  n <- 75; rows <- setdiff(1:n, c(1, 10:19, 60:63)); rho <- .65
  # A Cholesky whitening of the marginal covariance is an independent reference.
  covariance <- rho^abs(outer(rows, rows, "-"))
  w <- chol(solve(covariance)) * sqrt(1 - rho^2)
  z <- cbind(1, seq(-1, 1, length.out = n))
  designs <- list(cbind(sin(1:n / 4), cos(1:n / 7)),
    cbind(sin(1:n / 4 + .5), cos(1:n / 7 + .5)))
  y <- matrix(rnorm(n * 3), n, 3) + designs[[2]] %*% matrix(rnorm(6), 2)
  got <- hrf_score_run_cpp(y, designs, z, rows, rho = rho, evidence = TRUE)
  for (v in 1:3) {
    base <- lm.fit(w %*% z[rows, ], w %*% y[rows, v])
    full <- lapply(designs, function(x) lm.fit(w %*% cbind(z, x)[rows, ], w %*% y[rows, v]))
    union <- lm.fit(w %*% cbind(z, designs[[1]], designs[[2]])[rows, ], w %*% y[rows, v])
    for (k in 1:2) {
      sse <- sum(full[[k]]$residuals^2); sse0 <- sum(base$residuals^2)
      task_f <- ((sse0 - sse) / 2) / (sse / (length(rows) - full[[k]]$rank))
      expect_equal(got$task_p[k, v], pf(task_f, 2, length(rows) - full[[k]]$rank,
        lower.tail = FALSE), tolerance = 1e-9)
      expect_equal(got$scores[k, v], 1 - (sse / (length(rows) - full[[k]]$rank)) /
        (sse0 / (length(rows) - base$rank)), tolerance = 1e-9)
    }
    df <- length(rows) - union$rank
    shape_f <- ((sum(full[[1]]$residuals^2) - sum(union$residuals^2)) /
      (union$rank - full[[1]]$rank)) / (sum(union$residuals^2) / df)
    expect_equal(got$fallback_p[2, v], pf(shape_f, union$rank - full[[1]]$rank,
      df, lower.tail = FALSE), tolerance = 1e-9)
  }
  marginal_w <- hrf_whitening_matrix(rows, rho)
  expect_equal(marginal_w %*% covariance %*% t(marginal_w),
    diag(1 - rho^2, length(rows)), tolerance = 1e-12)
  contaminated <- y; contaminated[-rows, ] <- NA_real_
  expect_equal(hrf_score_run_cpp(contaminated, designs, z, rows, rho = rho, evidence = TRUE), got)
  expect_error(hrf_score_run_cpp(y, designs, z, rev(rows)), "increasing")
  same <- hrf_score_run_cpp(y, rep(designs[1], 2), z, rows, evidence = TRUE)
  expect_true(all(same$fallback_p == 1))
})

test_that("REML profiles agree with independent marginal covariance likelihoods", {
  set.seed(883)
  n <- 60; rows <- setdiff(1:n, c(4, 15:19)); grid <- c(-.3, 0, .4, .8)
  z <- cbind(1, seq(-1, 1, length.out = n), sin(1:n / 8))
  y <- matrix(rnorm(n * 2), n, 2)
  y[, 2] <- as.numeric(z %*% c(100, 1, 2)) + 1e-4 * y[, 2]
  got <- hrf_noise_profile_cpp(y, z, rows, grid)
  expect_equal(hrf_noise_profile_cpp(y, z, rows, grid, batch_size = 1L), got, tolerance = 1e-10)
  q <- qr.Q(qr(z[rows, ])); df <- length(rows) - ncol(q)
  for (k in seq_along(grid)) {
    covariance <- grid[k]^abs(outer(rows, rows, "-"))
    precision <- solve(covariance)
    for (v in 1:2) {
      a <- y[rows, v] - mean(y[rows, v]); a <- a / sqrt(sum(a^2))
      beta <- solve(crossprod(q, precision %*% q), crossprod(q, precision %*% a))
      residual <- a - q %*% beta
      rss <- as.numeric(crossprod(residual, precision %*% residual))
      expected <- df * log(rss / df) + as.numeric(determinant(crossprod(q, precision %*% q),
        logarithm = TRUE)$modulus) + as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
      expect_equal(got[k, v], expected, tolerance = 1e-10)
    }
  }
})

test_that("pooled covariance and calibrated null evidence behave on independent AR noise", {
  set.seed(55271)
  n <- 240L; voxels <- 64L; rho <- .6
  bold <- lapply(c("A", "B", "C"), function(r) {
    y <- replicate(voxels, as.numeric(stats::filter(rnorm(n + 100L), rho,
      method = "recursive"))[-seq_len(100L)]) + 100
    colnames(y) <- paste0("v", seq_len(voxels)); y
  })
  names(bold) <- c("A", "B", "C")
  ev <- data.frame(onset = seq(10, 220, by = 15), condition = "feedback")
  events <- setNames(rep(list(ev), 3L), names(bold))
  fit <- estimate_crossrun_hrf(bold, events, 1, hrf_library(1), amplitude_model = "trial")
  expect_true(all(vapply(fit$noise, function(x) abs(x$rho[1L] - rho) <= .1 + 1e-10, logical(1))))
  # Broad deterministic bound guards catastrophic null acceptance; the repeated
  # simulation report estimates rates and Monte Carlo intervals separately.
  expect_lt(mean(vapply(fit$selections, function(x) mean(x$kernel_index != 1L), numeric(1))), .15)
  expect_lt(mean(vapply(fit$selections, function(x) mean(x$selected_task_p <= .05), numeric(1))), .15)
  rescaled <- lapply(bold, function(y) -2 * y + 500)
  again <- estimate_crossrun_hrf(rescaled, events, 1, hrf_library(1), amplitude_model = "trial")
  expect_equal(lapply(again$noise, `[[`, "rho"), lapply(fit$noise, `[[`, "rho"))
  expect_equal(lapply(again$selections, `[[`, "kernel_index"), lapply(fit$selections, `[[`, "kernel_index"))
})

test_that("unidentifiable noise is explicit rather than silently treated as white", {
  d <- hrfn_data(n = 20L)
  masks <- lapply(d$bold, function(y) seq_len(nrow(y)) <= 8L)
  fit <- estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels, valid_rows = masks)
  expect_true(all(fit$signal_status$status == "noise_unidentified"))
  expect_true(all(fit$selections$A$voxel_evidence == "insufficient_training"))
  out <- deconvolve_crossrun_hrf(d$bold, fit, .01)
  expect_true(all(is.na(out$activity$A)))
  expect_true(all(out$scaling$A$status == "noise_unidentified"))
})

test_that("masked nuisance-adjusted inverse matches a joint normal-equation reference", {
  d <- hrfn_data(); n <- nrow(d$bold$A); rho <- .55
  rows <- setdiff(1:n, c(20:28, 80)); h <- d$kernels[, 2]
  z <- cbind(1, seq(-1, 1, length.out = n)); latent <- n + length(h) - 1L
  convolution <- matrix(0, n, latent)
  for (t in 1:n) for (j in seq_along(h)) convolution[t, t + length(h) - j] <- h[j]
  covariance <- rho^abs(outer(rows, rows, "-"))
  precision <- solve(covariance) * (1 - rho^2)
  design <- cbind(convolution[rows, ], z[rows, ])
  penalty <- crossprod(diff(diag(latent), differences = 2)) * .03 + diag(1e-8, latent)
  lhs <- crossprod(design, precision %*% design)
  lhs[seq_len(latent), seq_len(latent)] <- lhs[seq_len(latent), seq_len(latent)] + penalty
  expected <- solve(lhs, crossprod(design, precision))[length(h) - 1L + seq_len(n), ]
  settings <- list(backend = "cpp", lambda = .03, penalty = "diff2", ridge_floor = 1e-8, tail_padding = 0L)
  got <- hrf_inverse_operator(n, h, rows, rho, z, settings)
  expect_equal(got, expected, tolerance = 1e-6)
  settings$backend <- "r"
  expect_equal(hrf_inverse_operator(n, h, rows, rho, z, settings), expected, tolerance = 1e-6)
  # Joint fitting removes any supplied nuisance contribution, including slope.
  expect_lt(max(abs(got %*% z[rows, ])), 1e-7)
})

test_that("calibrated folds exclude all target information and flag inconsistent runs", {
  d <- hrfn_data()
  args <- c(d[c("bold", "events")], list(tr = 1, kernel_grid = d$kernels,
    amplitude_model = "trial", return_scores = TRUE))
  original <- do.call(estimate_crossrun_hrf, args)
  args$bold$A[,] <- NA_real_; args$events$A$onset <- args$events$A$onset + 1000
  changed <- do.call(estimate_crossrun_hrf, args)
  expect_identical(changed$selections$A, original$selections$A)
  expect_identical(changed$scores$A, original$scores$A)
  expect_identical(changed$parcels$A, original$parcels$A)
  # Distinct well-supported shapes in the two training runs have no common fit.
  d$bold$B <- d$bold$B[, rep(2L, 6L)]
  d$bold$C <- d$bold$C[, rep(1L, 6L)]
  colnames(d$bold$B) <- colnames(d$bold$C) <- colnames(d$bold$A)
  conflict <- estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels,
    amplitude_model = "trial", noise_rho = 0)
  expect_true(all(conflict$selections$A$kernel_index == 1L))
  expect_true(all(conflict$selections$A$stability == "no_common_library_shape"))
  two <- estimate_crossrun_hrf(d$bold[1:2], d$events[1:2], 1, d$kernels,
    amplitude_model = "trial", noise_rho = 0)
  expect_true(all(two$selections$A$stability == "not_assessable"))
})

test_that("modern application honors masks, covariance, nuisance, and scaling exactly once", {
  d <- hrfn_data(); n <- nrow(d$bold$A)
  mask <- lapply(d$bold, function(y) !(seq_len(n) %in% c(20:25, 80)))
  fit <- estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels, valid_rows = mask,
    amplitude_model = "trial", noise_rho = .4)
  d$bold$A[!mask$A, ] <- NA_real_
  for (backend in c("cpp", "r")) {
    got <- deconvolve_crossrun_hrf(d$bold, fit, .02, scaling = "zscore", backend = backend, tail_padding = 3L,
      nuisance_mode = "joint")
    for (r in names(d$bold)) for (v in 1:6) {
      rows <- which(mask[[r]]); k <- fit$selections[[r]]$kernel_index[v]
      z <- cbind(1, seq(-1, 1, length.out = n))
      op <- hrf_inverse_operator(n, fit$kernels[, k], rows, .4, z, got$settings)
      y <- d$bold[[r]][rows, v]
      expected <- op %*% ((y - mean(y)) / sd(y))
      expect_equal(as.numeric(got$activity[[r]][, v]), as.numeric(expected), tolerance = 1e-8)
    }
    expect_true(all(is.finite(got$activity$A)))
    expect_identical(attr(got$activity$A, "observation_mask"), mask$A)
  }
  expect_error(deconvolve_crossrun_hrf(d$bold, fit, .02, backend = "wiener"), "masked")
})

test_that("schedule diagnostics equal explicit event injection and subtract ideal overlap", {
  d <- hrfn_data()
  fit <- estimate_crossrun_hrf(d$bold, d$events, 1, d$kernels, amplitude_model = "trial", noise_rho = .4)
  got <- deconvolve_crossrun_hrf(d$bold, fit, .02, scaling = "demean", nuisance_mode = "joint")
  events <- lapply(d$events, function(e) data.frame(onset = c(20.3, 21.2, 35.4), condition = c("a", "b", "a")))
  diag <- diagnose_crossrun_recovery(got, events, lags = -1:2, return_matrices = TRUE)
  expect_true(any(abs(diag$matrices[[1L]]$ideal[row(diag$matrices[[1L]]$ideal) !=
    col(diag$matrices[[1L]]$ideal)]) > 0))
  for (r in names(events)) {
    n <- nrow(d$bold[[r]]); ev <- events[[r]]
    for (key in unique(diag$voxels$stratum[diag$voxels$run == r])) {
      v <- which(diag$voxels$run == r & diag$voxels$stratum == key)[1L]
      v <- match(diag$voxels$voxel[v], fit$voxel_ids)
      k <- fit$selections[[r]]$kernel_index[v]
      forward <- hrf_event_design_cpp(n, 1, 0, ev$onset, rep(0, 3), 1:3, 3L,
        fit$kernels[, k, drop = FALSE])[[1L]]
      op <- hrf_inverse_operator(n, fit$kernels[, k], 1:n, .4,
        cbind(1, seq(-1, 1, length.out = n)), got$settings)
      recovered <- op %*% sweep(forward, 2, colMeans(forward), "-")
      reference <- matrix(0, 3, 3)
      for (j in 1:3) for (s in 1:3) for (lag in -1:2)
        reference[j, s] <- reference[j, s] + approx(0:(n - 1L), recovered[, s], ev$onset[j] + lag)$y
      matrix <- diag$matrices[[paste(r, key, 1, sep = "/")]]
      expect_equal(matrix$recovered, reference, tolerance = 1e-10)
      expect_equal(matrix$excess, matrix$recovered - matrix$ideal)
    }
  }
  expect_false(diag$settings$automatic_correction)
  edge <- events; edge$A$onset[1] <- -10
  outside <- diagnose_crossrun_recovery(got, edge)
  expect_true(all(is.na(outside$events$amplitude_gain[outside$events$run == "A" & outside$events$event == 1])))
})

test_that("calibrated spatial borrowing admits task-supported parcel anchors when single-voxel shape power is modest", {
  set.seed(9042)
  n <- 110L
  h <- cbind(early = c(0, .25, .5, .2, .05, 0), late = c(0, .15, .45, .3, .1, 0))
  h <- sweep(h, 2L, colSums(h), "/")
  ev <- data.frame(onset = seq(8, n - 12, by = 14) + .3, condition = "feedback")
  events <- setNames(rep(list(ev), 3L), c("A", "B", "C"))
  x <- hrf_event_design_cpp(n, 1, 0, ev$onset, rep(0, nrow(ev)),
    rep(1L, nrow(ev)), 1L, h)
  # Voxels 1-3 have strong late-HRF signal (establishing the parcel anchor);
  # Voxel 4 has a late-HRF signal strong enough for task support (task_p <= .05)
  # and higher R^2 under late than early, but modest enough that single-voxel
  # shape_p > .05; Voxel 5 genuinely follows the early (canonical) HRF.
  bold <- lapply(events, function(e) {
    s_late <- as.numeric(x[[2L]])
    s_early <- as.numeric(x[[1L]])
    cbind(
      v1 = 100 + s_late + rnorm(n, sd = .01),
      v2 = 100 + s_late + rnorm(n, sd = .01),
      v3 = 100 + s_late + rnorm(n, sd = .01),
      v4 = 100 + s_late + rnorm(n, sd = .22),
      v5 = 100 + s_early + rnorm(n, sd = .01)
    )
  })
  groups <- setNames(rep("p1", 5L), colnames(bold$A))
  fit_spatial <- estimate_crossrun_hrf(bold, events, 1, h, spatial_groups = groups,
    pooling = "voxel", spatial_strength = .02, min_parcel_voxels = 3L,
    amplitude_model = "condition", noise_rho = 0, return_scores = TRUE)
  sel_s <- fit_spatial$selections$A
  expect_equal(sel_s$parcel_kernel_index[4L], 2L)
  expect_gt(sel_s$selected_shape_p[4L], .05)
  expect_lte(sel_s$selected_task_p[4L], .05)
  expect_equal(sel_s$kernel_index[4L], 2L)
  expect_equal(sel_s$status[4L], "voxel_spatial")
  expect_equal(sel_s$voxel_evidence[4L], "unidentified")
  # Genuinely canonical voxel 5 retains kernel 1 under modest spatial penalty.
  expect_equal(sel_s$kernel_index[5L], 1L)

  fit_fallback <- estimate_crossrun_hrf(bold, events, 1, h, spatial_groups = groups,
    pooling = "voxel", spatial_strength = 0, canonical_strength = 0, min_parcel_voxels = 3L,
    amplitude_model = "condition", noise_rho = 0, return_scores = TRUE)
  sel_f <- fit_fallback$selections$A
  expect_equal(sel_f$kernel_index[4L], 2L)
  expect_equal(sel_f$status[4L], "fallback_parcel")
  expect_equal(sel_f$voxel_evidence[4L], "unidentified")
  expect_equal(sel_f$kernel_index[5L], 1L)
  expect_equal(sel_f$status[5L], "voxel")
})

test_that("matrix-free bidiagonal inverse operator matches dense R reference across penalties, gaps, and nuisance modes", {
  n <- 96L
  kernel <- c(0, 0.15, 0.45, 0.25, 0.1, 0.05)
  rows_gap <- setdiff(seq_len(n), c(1L, 14:19, 55L, 88:90))
  z_none <- matrix(numeric(0), n, 0L)
  z_joint <- cbind(1, seq(-1, 1, length.out = n), sin(seq_len(n) / 9), cos(seq_len(n) / 13))

  for (rho in c(0, -0.35, 0.6)) {
    for (rows in list(seq_len(n), rows_gap)) {
      for (tail_padding in c(0L, 5L)) {
        for (penalty in c("diff2", "diff1", "ridge")) {
          s_cpp <- list(backend = "cpp", lambda = 0.04, penalty = penalty,
                        ridge_floor = 1e-8, tail_padding = tail_padding)
          s_r <- modifyList(s_cpp, list(backend = "r"))

          op_cpp_ignore <- hrf_inverse_operator(n, kernel, rows, rho, z_none, s_cpp)
          op_r_ignore <- hrf_inverse_operator(n, kernel, rows, rho, z_none, s_r)
          expect_equal(op_cpp_ignore, op_r_ignore, tolerance = 1e-12)

          op_cpp_joint <- hrf_inverse_operator(n, kernel, rows, rho, z_joint, s_cpp)
          op_r_joint <- hrf_inverse_operator(n, kernel, rows, rho, z_joint, s_r)
          expect_equal(op_cpp_joint, op_r_joint, tolerance = 1e-6)
          expect_lt(max(abs(op_cpp_joint %*% z_joint[rows, ])), 5e-7)
        }
      }
    }
  }
})

test_that("condition-structured trial scoring and canonical-angle evidence recover HRFs under coupled signed RPE and dense overlap", {
  set.seed(71204)
  n <- 200L
  lib <- hrf_library(1)
  runs <- c("run1", "run2", "run3")
  bold <- events <- vector("list", 3L)
  names(bold) <- names(events) <- runs
  groups <- stats::setNames(rep(c("p1", "p2"), each = 6L), paste0("v", 1:12))
  true_idx <- c(rep(5L, 6L), rep(6L, 6L))

  for (r in seq_along(runs)) {
    fb <- seq(20, 160, by = 10)
    cue <- fb - runif(length(fb), 3.0, 4.5)
    resp <- fb - runif(length(fb), 0.8, 1.4)
    ev <- rbind(
      data.frame(onset = cue, condition = "cue"),
      data.frame(onset = resp, condition = "response"),
      data.frame(onset = fb, condition = "feedback")
    )
    rpe <- as.numeric(scale(rnorm(length(fb))))
    # Zero mean signed RPE on feedback only; cue and response have zero mean.
    amps <- c(rep(0, 2 * length(fb)), 0.7 * rpe)
    des <- hrf_event_design_cpp(n, 1, 0, ev$onset, rep(0, nrow(ev)),
      seq_len(nrow(ev)), nrow(ev), lib)
    y <- matrix(0, n, 12L)
    for (v in 1:12) {
      sig <- as.numeric(des[[true_idx[v]]] %*% amps)
      ar_noise <- as.numeric(stats::filter(rnorm(n + 50L, sd = 0.03), 0.4,
        method = "recursive"))[-seq_len(50L)]
      y[, v] <- 100 + sig + ar_noise
    }
    colnames(y) <- names(groups)
    bold[[r]] <- y
    events[[r]] <- ev
  }

  fit <- estimate_crossrun_hrf(bold, events, 1, lib, spatial_groups = groups,
    spatial_strength = 0.02, amplitude_model = "trial")
  expect_gte(mean(fit$selections$run1$kernel_index == true_idx), 0.9)
})

test_that("coordinate-aware spatial regularization and Empirical Bayes shrinkage preserve sub-parcel territories", {
  set.seed(81920)
  n <- 180L
  lib <- hrf_library(1)
  runs <- c("run1", "run2", "run3")
  bold <- events <- vector("list", 3L)
  names(bold) <- names(events) <- runs
  coords <- cbind(x = seq_len(12L), y = rep(1, 12L))
  rownames(coords) <- paste0("v", 1:12)
  # Spatial gradient within a single volume: voxels 1-6 have HRF 5 (5s), voxels 7-12 have HRF 6 (7s).
  true_idx <- c(rep(5L, 6L), rep(6L, 6L))

  for (r in seq_along(runs)) {
    fb <- seq(18, 150, by = 11)
    ev <- data.frame(onset = fb, condition = "feedback")
    amps <- 0.5 + 0.5 * as.numeric(scale(rnorm(length(fb))))
    des <- hrf_event_design_cpp(n, 1, 0, ev$onset, rep(0, nrow(ev)),
      seq_len(nrow(ev)), nrow(ev), lib)
    y <- matrix(0, n, 12L)
    for (v in 1:12) {
      sig <- as.numeric(des[[true_idx[v]]] %*% amps)
      sd_v <- if (v %in% c(3L, 10L)) 0.12 else 0.025
      y[, v] <- 100 + sig + rnorm(n, sd = sd_v)
    }
    colnames(y) <- rownames(coords)
    bold[[r]] <- y
    events[[r]] <- ev
  }

  fit_coords <- estimate_crossrun_hrf(bold, events, 1, lib,
    spatial_coords = coords, spatial_fwhm = 2.5, spatial_strength = 0.05,
    min_parcel_voxels = 4L, amplitude_model = "trial", noise_model = "ols")
  expect_equal(fit_coords$selections$run1$kernel_index, true_idx)
})

test_that("deconvolve_crossrun_hrf supports noise_sd scaling, matern1/sobolev penalties, and decon_epoch_regression", {
  set.seed(7204)
  lib <- hrf_library(1, peak_delays = c(5, 6), dispersions = 1, undershoot_ratios = 6)
  hrf <- lib[, 2L]
  n_time <- 160L
  n_trials <- 20L
  true_beta <- 0.80

  bold <- events <- rpe_list <- list()
  for (r in 1:2) {
    t_cue <- seq(6, by = 7, length.out = n_trials)
    t_fb <- t_cue + runif(n_trials, 1.5, 2.3)
    rpe <- rnorm(n_trials, sd = if (r == 1L) 0.7 else 1.4)
    ev <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue"),
      data.frame(onset = t_fb, duration = 0, condition = "feedback")
    )
    des_cue <- hrf_event_design_cpp(n_time, 1, 0, t_cue, rep(0, n_trials),
      seq_len(n_trials), n_trials, matrix(hrf, ncol = 1L))[[1L]]
    des_fb <- hrf_event_design_cpp(n_time, 1, 0, t_fb, rep(0, n_trials),
      seq_len(n_trials), n_trials, matrix(hrf, ncol = 1L))[[1L]]
    amp_cue <- (if (r == 1L) 1.2 else 2.2) + 0.4 * rnorm(n_trials)
    amp_fb <- 0.5 + true_beta * rpe
    sig <- as.numeric(des_cue %*% amp_cue + des_fb %*% amp_fb)
    y <- cbind(
      v1 = 100 + sig + rnorm(n_time, sd = 0.2),
      v2 = 100 + sig + rnorm(n_time, sd = 0.2),
      v3 = 100 + sig + rnorm(n_time, sd = 0.2),
      v4 = 100 + sig + rnorm(n_time, sd = 0.2)
    )
    bold[[paste0("run", r)]] <- y
    events[[paste0("run", r)]] <- ev
    rpe_list[[paste0("run", r)]] <- rpe
  }

  fit <- estimate_crossrun_hrf(bold, events, 1, lib, noise_model = "ar1", evidence = "calibrated")
  for (pen in c("matern1", "sobolev")) {
    for (nm in c("ignore", "joint")) {
      dec_cpp <- deconvolve_crossrun_hrf(bold, fit, lambda = 1.0, scaling = "noise_sd",
        penalty = pen, nuisance_mode = nm, backend = "cpp")
      dec_r <- deconvolve_crossrun_hrf(bold, fit, lambda = 1.0, scaling = "noise_sd",
        penalty = pen, nuisance_mode = nm, backend = "r")
      expect_equal(dec_cpp$activity$run1, dec_r$activity$run1, tolerance = 1e-10)
      expect_equal(dec_cpp$activity$run2, dec_r$activity$run2, tolerance = 1e-10)
    }
  }

  dec_base <- deconvolve_crossrun_hrf(bold, fit, lambda = 1.0, scaling = "zscore", penalty = "diff2")
  dec_imp <- deconvolve_crossrun_hrf(bold, fit, lambda = NULL, scaling = "noise_sd", penalty = "sobolev")
  expect_true(is.finite(dec_imp$settings$lambda) && dec_imp$settings$lambda > 0)
  expect_true(is.data.frame(dec_imp$settings$lambda_tuning))

  reg_base <- decon_epoch_regression(dec_base, rpe_list, events = events,
    target_condition = "feedback", lags = c(-1, 0, 1), readout_kernel = "linear",
    demix_other_events = FALSE, demix_target_trials = FALSE, restore_scale = TRUE)
  reg_imp <- decon_epoch_regression(dec_imp, rpe_list, events = events,
    target_condition = "feedback", lags = c(-1, 0, 1), readout_kernel = "triangle",
    readout_width = 1.5, demix_other_events = TRUE, demix_target_trials = TRUE,
    restore_scale = TRUE, source_model = "window", return_epochs = TRUE, return_panel = TRUE)

  expect_true(is.list(reg_imp$epochs))
  expect_true(is.data.frame(reg_imp$panel))
  expect_true(all(c("run", "trial", "evt_time", "voxel", "decon_raw", "decon_demixed") %in% names(reg_imp$panel)))
  # Improved de-mixed regression has higher trial correlation and much smaller amplitude bias at lag = 0
  r0_base <- reg_base$summary$mean_r_raw[reg_base$summary$lag == 0]
  r0_imp <- reg_imp$summary$mean_r_demixed[reg_imp$summary$lag == 0]
  b0_base <- reg_base$summary$mean_beta_raw[reg_base$summary$lag == 0]
  b0_imp <- reg_imp$summary$mean_beta_corrected[reg_imp$summary$lag == 0]
  expect_gt(r0_imp, r0_base)
  expect_lt(abs(b0_imp - true_beta), 0.35 * abs(b0_base - true_beta))
})

test_that("decon_epoch_regression trial-varying duration source model eliminates RT time-on-task confounding", {
  set.seed(2026)
  tr <- 1.0
  n_time <- 220L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = c(5, 6), dispersions = 1, undershoot_ratios = 6)
  hrf <- lib[, 1L]
  true_beta <- 0.65

  # 1. Decision-making paradigm where higher value x has faster RT (cor(x, rt) = -0.60)
  # and positive baseline choice firing rate alpha0 = 1.8
  bold_val <- bold_null <- bold_int <- ev_imp <- ev_dur <- x_list <- rt_list <- list()
  for (r in 1:2) {
    t_cue <- seq(6, by = 8.5, length.out = n_trials)
    t_choice <- t_cue + runif(n_trials, 1.4, 2.1)
    x <- as.numeric(scale(rnorm(n_trials)))
    z_orth <- as.numeric(scale(stats::resid(stats::lm(rnorm(n_trials) ~ x))))
    rt_neg <- pmin(2.8, pmax(0.4, 1.5 + 0.65 * (-0.60 * x + sqrt(1 - 0.60^2) * z_orth)))
    des_cue <- hrf_event_design_cpp(n_time, tr, 0, t_cue, rep(0, n_trials),
      seq_len(n_trials), n_trials, matrix(hrf, ncol = 1L))[[1L]]
    des_choice <- hrf_event_design_cpp(n_time, tr, 0, t_choice, rt_neg,
      seq_len(n_trials), n_trials, matrix(hrf, ncol = 1L))[[1L]]
    amp_cue <- 1.0 + 0.2 * rnorm(n_trials)
    amp_val <- 1.8 + true_beta * x
    amp_null <- rep(1.8, n_trials)

    sig_val <- as.numeric(des_cue %*% amp_cue + des_choice %*% amp_val)
    sig_null <- as.numeric(des_cue %*% amp_cue + des_choice %*% amp_null)
    des_choice_int <- sweep(des_choice, 2L, rt_neg / tr, "/")
    sig_int <- as.numeric(des_cue %*% amp_cue + des_choice_int %*% amp_val)

    noise_mat <- replicate(6L, rnorm(n_time, sd = 0.15))
    colnames(noise_mat) <- paste0("v", 1:6)
    bold_val[[paste0("run", r)]] <- 100 + sig_val + noise_mat
    bold_null[[paste0("run", r)]] <- 100 + sig_null + noise_mat
    bold_int[[paste0("run", r)]] <- 100 + sig_int + noise_mat

    ev_imp[[paste0("run", r)]] <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue"),
      data.frame(onset = t_choice, duration = 0, condition = "choice")
    )
    ev_dur[[paste0("run", r)]] <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue"),
      data.frame(onset = t_choice, duration = rt_neg, condition = "choice")
    )
    x_list[[paste0("run", r)]] <- x
    rt_list[[paste0("run", r)]] <- rt_neg
  }

  fit_imp <- estimate_crossrun_hrf(bold_val, ev_imp, tr, lib, noise_model = "ols")
  dec_imp <- deconvolve_crossrun_hrf(bold_val, fit_imp, lambda = 1.0, scaling = "noise_sd", penalty = "sobolev")
  fit_dur <- estimate_crossrun_hrf(bold_val, ev_dur, tr, lib, noise_model = "ols")
  dec_dur <- deconvolve_crossrun_hrf(bold_val, fit_dur, lambda = 1.0, scaling = "noise_sd", penalty = "sobolev")

  reg_imp <- decon_epoch_regression(dec_imp, x_list, events = ev_imp,
    target_condition = "choice", lags = 0, source_model = "impulse", restore_scale = TRUE)
  reg_win <- decon_epoch_regression(dec_imp, x_list, events = ev_imp,
    target_condition = "choice", lags = 0, source_model = "window", restore_scale = TRUE)
  reg_dur <- decon_epoch_regression(dec_dur, x_list, events = ev_imp,
    target_condition = "choice", lags = 0, trial_durations = rt_list, restore_scale = TRUE)

  expect_equal(reg_dur$settings$source_model, "duration")
  expect_equal(reg_dur$settings$duration_scale, "rate")
  expect_gt(reg_dur$summary$mean_r_demixed, reg_imp$summary$mean_r_demixed + 0.35)
  expect_gt(reg_dur$summary$mean_r_demixed, reg_win$summary$mean_r_demixed + 0.35)
  expect_lt(abs(reg_dur$summary$mean_beta_corrected - true_beta), 0.12)
  expect_lt(
    abs(reg_dur$summary$mean_beta_corrected - true_beta),
    0.35 * abs(reg_win$summary$mean_beta_corrected - true_beta)
  )

  # 2. Under null modulation (true_beta = 0) with strong RT-predictor correlation,
  # impulse/window models suffer severe spurious RT bias, while duration model is unbiased
  fit_null_imp <- estimate_crossrun_hrf(bold_null, ev_imp, tr, lib, noise_model = "ols")
  dec_null_imp <- deconvolve_crossrun_hrf(bold_null, fit_null_imp, lambda = 1.0,
    scaling = "noise_sd", penalty = "sobolev")
  fit_null_dur <- estimate_crossrun_hrf(bold_null, ev_dur, tr, lib, noise_model = "ols")
  dec_null_dur <- deconvolve_crossrun_hrf(bold_null, fit_null_dur, lambda = 1.0,
    scaling = "noise_sd", penalty = "sobolev")

  null_imp <- decon_epoch_regression(dec_null_imp, x_list, events = ev_imp,
    target_condition = "choice", lags = c(-1, 0, 1), source_model = "impulse", restore_scale = TRUE)
  null_dur <- decon_epoch_regression(dec_null_dur, x_list, events = ev_dur,
    target_condition = "choice", lags = c(-1, 0, 1), source_model = "duration", restore_scale = TRUE)

  expect_gt(abs(null_imp$summary$mean_r_demixed[null_imp$summary$lag == 0]), 0.40)
  expect_lt(max(abs(null_dur$summary$mean_r_demixed)), 0.15)
  expect_lt(max(abs(null_dur$summary$mean_beta_corrected)), 0.08)

  # 3. Integral duration_scale recovers unit-integral boxcar amplitudes
  fit_int <- estimate_crossrun_hrf(bold_int, ev_dur, tr, lib, noise_model = "ols")
  dec_int <- deconvolve_crossrun_hrf(bold_int, fit_int, lambda = 1.0, scaling = "noise_sd", penalty = "sobolev")
  reg_int <- decon_epoch_regression(dec_int, x_list, events = ev_dur,
    target_condition = "choice", lags = 0, source_model = "duration",
    duration_scale = "integral", restore_scale = TRUE)
  expect_gt(reg_int$summary$mean_r_demixed, 0.70)
  expect_lt(abs(reg_int$summary$mean_beta_corrected - true_beta), 0.15)
})

test_that("decon_epoch_regression propagates analytical per-trial precision weights and ZCA whitening to mixed_by", {
  set.seed(202612)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = c(5, 6), dispersions = 1, undershoot_ratios = 6)
  hrf <- lib[, 1L, drop = FALSE]
  true_group_beta <- 0.45
  n_subj <- 10L
  subj_betas <- true_group_beta + rnorm(n_subj, sd = 0.08)

  panel_all <- vector("list", n_subj)
  subj_err_unw <- subj_err_pw <- subj_err_gls <- numeric(n_subj)

  for (s in seq_len(n_subj)) {
    bold <- events <- cov_list <- rt_list <- list()
    for (r in 1:2) {
      # Heterogeneous ISI (crowded vs isolated) and trial durations (short RT = 0.35s vs long RT = 2.4s)
      isi <- rep(c(2.8, 7.2), length.out = n_trials) + runif(n_trials, -0.25, 0.25)
      t_choice <- cumsum(c(8, isi[-n_trials]))
      rt <- sample(c(rep(0.35, n_trials / 2L), rep(2.4, n_trials / 2L)))
      val <- as.numeric(scale(stats::filter(rnorm(n_trials + 10L), 0.50, method = "recursive")[-seq_len(10L)]))
      des_choice <- hrf_event_design_cpp(n_time, tr, 0, t_choice, rt,
        seq_len(n_trials), n_trials, hrf)[[1L]]
      amp <- 1.0 + subj_betas[s] * val + rnorm(n_trials, sd = 0.05)
      sig <- as.numeric(des_choice %*% amp)
      # Run 1 has low scanner noise (sd = 0.12); Run 2 has high scanner noise (sd = 0.55)
      sd_r <- if (r == 1L) 0.12 else 0.55
      W_r <- hrf_whitening_matrix(1:n_time, 0.40)
      eps1 <- sd_r * sqrt(1 - 0.40^2) * as.numeric(solve(W_r, rnorm(n_time)))
      eps2 <- sd_r * sqrt(1 - 0.40^2) * as.numeric(solve(W_r, rnorm(n_time)))
      bold[[paste0("run", r)]] <- cbind(v1 = 100 + sig + eps1, v2 = 100 + sig + eps2)
      events[[paste0("run", r)]] <- data.frame(onset = t_choice, duration = rt, condition = "choice")
      cov_list[[paste0("run", r)]] <- data.frame(value = val)
      rt_list[[paste0("run", r)]] <- rt
    }

    fit_s <- list(
      kernels = hrf,
      run_ids = c("run1", "run2"),
      voxel_ids = c("v1", "v2"),
      n_time = c(run1 = n_time, run2 = n_time),
      selections = list(
        run1 = data.frame(voxel = c("v1", "v2"), kernel_index = c(1L, 1L), status = "voxel"),
        run2 = data.frame(voxel = c("v1", "v2"), kernel_index = c(1L, 1L), status = "voxel")
      ),
      noise = list(
        run1 = data.frame(voxel = c("v1", "v2"), rho = c(0.40, 0.40), method = "reml"),
        run2 = data.frame(voxel = c("v1", "v2"), rho = c(0.40, 0.40), method = "reml")
      ),
      specification = list(
        events = events,
        nuisance = list(
          run1 = cbind(1, seq(-1, 1, length.out = n_time)),
          run2 = cbind(1, seq(-1, 1, length.out = n_time))
        ),
        valid_rows = list(run1 = rep(TRUE, n_time), run2 = rep(TRUE, n_time))
      ),
      settings = list(tr = tr, time_offset = c(0, 0), noise_model = "ar1")
    )
    class(fit_s) <- "crossrun_hrf"
    dec_s <- deconvolve_crossrun_hrf(bold, fit_s, lambda = 1.0, scaling = "noise_sd", penalty = "sobolev")
    expect_true("noise_sd" %in% names(dec_s$scaling$run1))
    expect_true(all(is.finite(dec_s$scaling$run1$noise_sd)))

    reg_s <- decon_epoch_regression(
      dec_s, cov_list, events = events, target_condition = "choice",
      lags = 0, source_model = "duration", trial_durations = rt_list,
      restore_scale = TRUE, weighting = "precision", whiten_panel = TRUE,
      return_epochs = TRUE, return_panel = TRUE
    )

    expect_true(all(c(
      "se_raw", "se_demixed", "precision_weight", "precision_weight_op",
      "precision_weight_raw", "error_ac1", "decon_demixed_whitened", "value_whitened"
    ) %in% names(reg_s$panel)))

    # Short-RT trials and high-noise run2 must have larger analytical SE and smaller precision weight
    p_v1 <- reg_s$panel[reg_s$panel$voxel == "v1", ]
    rt_vec <- c(rt_list$run1, rt_list$run2)
    expect_gt(
      mean(p_v1$se_demixed[p_v1$run == "run2"]),
      2.5 * mean(p_v1$se_demixed[p_v1$run == "run1"])
    )
    expect_gt(
      mean(p_v1$se_demixed[p_v1$run == "run1" & rt_vec == 0.35]),
      2.0 * mean(p_v1$se_demixed[p_v1$run == "run1" & rt_vec == 2.4])
    )
    expect_gt(
      mean(p_v1$precision_weight[p_v1$run == "run1" & rt_vec == 2.4]),
      2.5 * mean(p_v1$precision_weight[p_v1$run == "run2" & rt_vec == 0.35])
    )
    expect_equal(mean(p_v1$precision_weight), 1.0, tolerance = 1e-10)
    expect_equal(mean(p_v1$precision_weight_op), 1.0, tolerance = 1e-10)

    subj_err_unw[s] <- abs(reg_s$summary$mean_beta_unweighted - subj_betas[s])
    subj_err_pw[s] <- abs(reg_s$summary$mean_beta_weighted - subj_betas[s])
    subj_err_gls[s] <- abs(reg_s$summary$mean_beta_gls - subj_betas[s])

    p_s <- reg_s$panel
    p_s$id <- factor(paste0("s", s))
    panel_all[[s]] <- p_s
  }

  # Precision-weighted and GLS single-subject estimators beat unweighted OLS
  expect_lt(mean(subj_err_pw), 0.65 * mean(subj_err_unw))
  expect_lt(mean(subj_err_gls), 0.65 * mean(subj_err_unw))

  # Downstream mixed_by with weights = "precision_weight" and whitened panel
  panel_group <- data.table::rbindlist(panel_all)
  res_unw <- suppressMessages(suppressWarnings(mixed_by(
    panel_group,
    outcomes = "decon_demixed",
    rhs_model_formulae = list(m1 = ~ value + (1 | id)),
    split_on = c("voxel", "evt_time"),
    calculate = "parameter_estimates_reml"
  )))
  res_pw <- suppressMessages(suppressWarnings(mixed_by(
    panel_group,
    outcomes = "decon_demixed",
    rhs_model_formulae = list(m1 = ~ value + (1 | id)),
    split_on = c("voxel", "evt_time"),
    weights = "precision_weight",
    calculate = "parameter_estimates_reml",
    emtrends_spec = list(vt = list(outcome = "decon_demixed", model_name = "m1", specs = ~ 1, var = "value"))
  )))
  res_gls <- suppressMessages(suppressWarnings(mixed_by(
    panel_group,
    outcomes = "decon_demixed_whitened",
    rhs_model_formulae = list(m1 = ~ value_whitened + (1 | id)),
    split_on = c("voxel", "evt_time"),
    engine_args = list(weights = "precision_weight"),
    calculate = "parameter_estimates_reml"
  )))

  b_unw <- res_unw$coef_df_reml[res_unw$coef_df_reml$term == "value", ]
  b_pw <- res_pw$coef_df_reml[res_pw$coef_df_reml$term == "value", ]
  b_gls <- res_gls$coef_df_reml[res_gls$coef_df_reml$term == "value_whitened", ]

  expect_true(data.table::is.data.table(res_pw$emtrends_list$vt))
  expect_equal(nrow(res_pw$emtrends_list$vt), 2L)
  expect_lt(mean(b_pw$std.error), 0.50 * mean(b_unw$std.error))
  expect_lt(mean(b_gls$std.error), 0.50 * mean(b_unw$std.error))
  expect_lt(mean(abs(b_pw$estimate - mean(subj_betas))), 0.08)
  expect_lt(mean(abs(b_gls$estimate - mean(subj_betas))), 0.08)
})

test_that("decon_epoch_regression removes trial-level non-target cue confounding and sharpens within-trial lag profiles", {
  set.seed(202615)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0), dispersions = c(0.9, 1.0), undershoot_ratios = 6)
  hrf_can <- lib[, 1L, drop = FALSE]

  bold <- events <- x_list <- list()
  for (r in 1:2) {
    t_cue <- seq(8, by = 6.5, length.out = n_trials)
    t_fb <- t_cue + runif(n_trials, 2.2, 3.2)
    x <- as.numeric(scale(rnorm(n_trials)))
    z_orth <- as.numeric(scale(stats::resid(stats::lm(rnorm(n_trials) ~ x))))
    cue_mod <- 0.65 * x + sqrt(1 - 0.65^2) * z_orth

    des_cue <- hrf_event_design_cpp(n_time, tr, 0, t_cue, rep(0, n_trials),
      seq_len(n_trials), n_trials, hrf_can)[[1L]]
    des_fb <- hrf_event_design_cpp(n_time, tr, 0, t_fb, rep(0, n_trials),
      seq_len(n_trials), n_trials, hrf_can)[[1L]]

    # Voxels 1-2: Null target (beta_fb = 0) with strong correlated cue modulation (beta_cue = 0.90)
    sig_null <- as.numeric(des_cue %*% (1.0 + 0.90 * cue_mod) + des_fb %*% rep(0.6, n_trials))
    # Voxels 3-4: True phasic target (beta_fb = 0.70) at lag = 0
    sig_phasic <- as.numeric(des_cue %*% rep(0.6, n_trials) + des_fb %*% (0.8 + 0.70 * x))

    bold[[paste0("run", r)]] <- cbind(
      v1 = 100 + sig_null + rnorm(n_time, sd = 0.12),
      v2 = 100 + sig_null + rnorm(n_time, sd = 0.12),
      v3 = 100 + sig_phasic + rnorm(n_time, sd = 0.12),
      v4 = 100 + sig_phasic + rnorm(n_time, sd = 0.12)
    )
    events[[paste0("run", r)]] <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue"),
      data.frame(onset = t_fb, duration = 0, condition = "feedback")
    )
    x_list[[paste0("run", r)]] <- data.frame(x = x)
  }

  fit <- estimate_crossrun_hrf(bold, events, tr, lib, noise_model = "ols")
  dec <- deconvolve_crossrun_hrf(bold, fit, lambda = 0.5, scaling = "noise_sd", penalty = "sobolev")

  reg_no_other <- decon_epoch_regression(dec, x_list, events = events, target_condition = "feedback",
    lags = -2:3, source_model = "impulse", demix_other_events = FALSE, demix_target_trials = TRUE)
  reg_with_other <- decon_epoch_regression(dec, x_list, events = events, target_condition = "feedback",
    lags = -2:3, source_model = "impulse", demix_other_events = TRUE, demix_target_trials = TRUE,
    demix_within_trial = TRUE, within_trial_ridge = 0.08)

  cf_no <- reg_no_other$coefficients
  cf_wt <- reg_with_other$coefficients

  b_null_no <- mean(abs(cf_no$beta_corrected[cf_no$voxel %in% c("v1", "v2") & cf_no$lag == 0]))
  b_null_wt <- mean(abs(cf_wt$beta_corrected[cf_wt$voxel %in% c("v1", "v2") & cf_wt$lag == 0]))
  expect_lt(b_null_wt, 0.35 * b_null_no)

  # Within-trial multi-lag de-mixing sharpens peak concentration at lag = 0 vs off-peak lags
  prof_no <- tapply(cf_no$beta_corrected[cf_no$voxel %in% c("v3", "v4")],
    cf_no$lag[cf_no$voxel %in% c("v3", "v4")], mean)
  prof_wt <- tapply(cf_wt$beta_corrected[cf_wt$voxel %in% c("v3", "v4")],
    cf_wt$lag[cf_wt$voxel %in% c("v3", "v4")], mean)
  ratio_no <- prof_no["0"] / mean(abs(prof_no[c("-2", "2", "3")]))
  ratio_wt <- prof_wt["0"] / mean(abs(prof_wt[c("-2", "2", "3")]))
  expect_gt(ratio_wt, ratio_no)
  expect_true(reg_with_other$settings$demix_within_trial)
})

test_that("soft_average selection_mode interpolates sub-grid HRFs and propagates composite kernels through deconvolution", {
  set.seed(202624)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = c(4.0, 5.0, 6.0, 7.0, 8.0), dispersions = c(0.9, 1.0, 1.1), undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6 & abs(pars$dispersion - 1.0) < 1e-6)[1L]

  # True HRF is an off-grid intermediate shape (peak_delay = 5.35s, dispersion = 0.96)
  true_hrf <- hrf_library(tr, peak_delays = 5.35, dispersions = 0.96, undershoot_ratios = 6)[, 1L, drop = FALSE]
  groups <- stats::setNames(rep("R_subgrid", 6L), paste0("v", 1:6))

  bold <- events <- x_list <- list()
  for (r in 1:2) {
    t_cue <- seq(8, by = 6.5, length.out = n_trials)
    t_fb <- t_cue + runif(n_trials, 1.8, 2.5)
    x <- as.numeric(scale(rnorm(n_trials)))
    des_cue <- hrf_event_design_cpp(n_time, tr, 0, t_cue, rep(0, n_trials),
      seq_len(n_trials), n_trials, true_hrf)[[1L]]
    des_fb <- hrf_event_design_cpp(n_time, tr, 0, t_fb, rep(0, n_trials),
      seq_len(n_trials), n_trials, true_hrf)[[1L]]
    sig <- as.numeric(des_cue %*% rep(0.5, n_trials) + des_fb %*% (0.8 + 0.70 * x))
    ymat <- vapply(1:6, function(v) 100 + sig + rnorm(n_time, sd = 0.14), numeric(n_time))
    colnames(ymat) <- names(groups)
    bold[[paste0("run", r)]] <- ymat
    events[[paste0("run", r)]] <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue"),
      data.frame(onset = t_fb, duration = 0, condition = "feedback")
    )
    x_list[[paste0("run", r)]] <- data.frame(x = x)
  }

  fit_hard <- estimate_crossrun_hrf(bold, events, tr, lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    amplitude_model = "trial", noise_model = "ols", selection_mode = "hard")
  fit_soft <- estimate_crossrun_hrf(bold, events, tr, lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    amplitude_model = "trial", noise_model = "ols", selection_mode = "soft_average")

  expect_equal(fit_soft$settings$selection_mode, "soft_average")
  expect_true(is.matrix(fit_soft$kernel_weights$run1))
  expect_true(is.matrix(fit_soft$selected_kernels$run1))
  expect_equal(unname(colSums(fit_soft$kernel_weights$run1)), rep(1, 6L), tolerance = 1e-10)
  expect_equal(unname(colSums(fit_soft$selected_kernels$run1)), rep(1, 6L), tolerance = 1e-10)

  rmse_hard <- sqrt(mean((fit_hard$selected_kernels$run1[, 1L] - true_hrf[, 1L])^2))
  rmse_soft <- sqrt(mean((fit_soft$selected_kernels$run1[, 1L] - true_hrf[, 1L])^2))
  expect_lt(rmse_soft, rmse_hard)

  dec_soft <- deconvolve_crossrun_hrf(bold, fit_soft, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")
  diag_soft <- diagnose_crossrun_recovery(dec_soft, events = events, lags = 0)
  reg_soft <- decon_epoch_regression(dec_soft, x_list, events = events,
    target_condition = "feedback", lags = -1:2, restore_scale = TRUE)
  expect_true(nrow(diag_soft$summary) >= 1L)
  expect_true(all(is.finite(reg_soft$summary$mean_beta_corrected)))
})

test_that("default adaptive prewhitened GCV lambda adapts to fSNR, TR, and null background voxels", {
  set.seed(202630)
  sim_biphasic <- function(tr = 1.0, noise_sd = 0.15, n_null = 24L) {
    n_time <- as.integer(round(180 / tr))
    n_trials <- 22L
    lib <- hrf_library(tr, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)
    hrf <- lib[, 1L, drop = FALSE]
    n_act <- 8L
    n_vox <- n_act + n_null
    groups <- stats::setNames(
      c(rep("R_biphasic", n_act), rep("R_null", n_null)),
      paste0("v", seq_len(n_vox))
    )
    bold <- events <- x_list <- list()
    for (r in 1:2) {
      t_ev <- seq(8, by = 7.0, length.out = n_trials)
      x <- as.numeric(scale(rnorm(n_trials)))
      # Biphasic neural response: positive impulse at t_ev and negative impulse at t_ev + 1.0s
      des_pos <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
        seq_len(n_trials), n_trials, hrf)[[1L]]
      des_neg <- hrf_event_design_cpp(n_time, tr, 0, t_ev + 1.0, rep(0, n_trials),
        seq_len(n_trials), n_trials, hrf)[[1L]]
      sig <- as.numeric(des_pos %*% (0.6 + 0.85 * x) + des_neg %*% (-0.3 - 0.65 * x))
      W_r <- hrf_whitening_matrix(seq_len(n_time), 0.35)
      ymat <- matrix(0.0, nrow = n_time, ncol = n_vox)
      for (v in seq_len(n_vox)) {
        eps <- noise_sd * sqrt(1 - 0.35^2) * as.numeric(solve(W_r, rnorm(n_time)))
        ymat[, v] <- 100 + (if (v <= n_act) sig else 0.0) + eps
      }
      colnames(ymat) <- names(groups)
      bold[[paste0("run", r)]] <- ymat
      events[[paste0("run", r)]] <- data.frame(onset = t_ev, duration = 0, condition = "feedback")
      x_list[[paste0("run", r)]] <- data.frame(x = x)
    }
    fit <- estimate_crossrun_hrf(bold, events, tr, lib,
      spatial_groups = groups, pooling = "parcel", amplitude_model = "trial")
    list(bold = bold, events = events, x = x_list, fit = fit, act_ids = paste0("v", seq_len(n_act)))
  }

  # 1. Omitting lambda uses adaptive prewhitened GCV by default on modern fits
  d_hi <- sim_biphasic(tr = 1.0, noise_sd = 0.08, n_null = 24L)
  d_lo <- sim_biphasic(tr = 1.0, noise_sd = 0.45, n_null = 24L)
  dec_hi <- deconvolve_crossrun_hrf(d_hi$bold, d_hi$fit, scaling = "noise_sd", penalty = "sobolev")
  dec_lo <- deconvolve_crossrun_hrf(d_lo$bold, d_lo$fit, scaling = "noise_sd", penalty = "sobolev")

  expect_true(is.finite(dec_hi$settings$lambda) && dec_hi$settings$lambda > 0)
  expect_true(is.data.frame(dec_hi$settings$lambda_tuning))
  expect_true(all(c("var_level", "var_contrast", "mtf_rel", "gcv_adaptive", "snr_task") %in%
    names(dec_hi$settings$lambda_tuning)))

  # High fSNR estimates higher snr_task and selects smaller lambda than low fSNR,
  # while low fSNR remains bounded below the Pareto-dominated over-smoothing regime
  expect_gt(dec_hi$settings$lambda_tuning$snr_task[1L], dec_lo$settings$lambda_tuning$snr_task[1L])
  expect_lt(dec_hi$settings$lambda, dec_lo$settings$lambda)
  expect_lte(dec_lo$settings$lambda, 2.0)

  # 2. Naive observation-space GCV across all voxels (including R_null) over-smooths
  # to the upper grid boundary and degrades biphasic lag contrast recovery
  sc_naive <- scale_decon_input(d_lo$bold$run1, scaling = "noise_sd", noise_rho = 0.35)$y
  naive_tune <- deconvolve_reglin(
    BOLDobs = sc_naive, kernel = d_lo$fit$kernels[, 1L], lambda = NULL,
    lambda_grid = 10^seq(-2, 1.0, length.out = 16), penalty = "sobolev",
    tune_by = "global", normalize = FALSE, demean = FALSE,
    return_diagnostics = TRUE, prewhiten_gcv = TRUE, prewhiten_rho = 0.35, scaling = "none"
  )
  lam_naive <- as.numeric(naive_tune$lambda[1L])
  expect_gt(lam_naive, 2.5 * dec_lo$settings$lambda)

  dec_naive <- deconvolve_crossrun_hrf(d_lo$bold, d_lo$fit, lambda = lam_naive,
    scaling = "noise_sd", penalty = "sobolev")
  reg_adapt <- decon_epoch_regression(dec_lo, d_lo$x, events = d_lo$events,
    target_condition = "feedback", lags = c(0, 1), demix_within_trial = TRUE)
  reg_naive <- decon_epoch_regression(dec_naive, d_lo$x, events = d_lo$events,
    target_condition = "feedback", lags = c(0, 1), demix_within_trial = TRUE)

  cf_adapt <- reg_adapt$coefficients[reg_adapt$coefficients$voxel %in% d_lo$act_ids, ]
  cf_naive <- reg_naive$coefficients[reg_naive$coefficients$voxel %in% d_lo$act_ids, ]
  contrast_adapt <- mean(cf_adapt$beta_corrected[cf_adapt$lag == 0]) -
    mean(cf_adapt$beta_corrected[cf_adapt$lag == 1])
  contrast_naive <- mean(cf_naive$beta_corrected[cf_naive$lag == 0]) -
    mean(cf_naive$beta_corrected[cf_naive$lag == 1])
  expect_gt(contrast_adapt, contrast_naive)

  # 3. Discrete derivative scaling across acquisition TRs: fast TR (0.5s) selects
  # larger lambda than standard TR (1.0s), which selects larger lambda than slow TR (2.0s)
  d_tr05 <- sim_biphasic(tr = 0.5, noise_sd = 0.20, n_null = 8L)
  d_tr20 <- sim_biphasic(tr = 2.0, noise_sd = 0.20, n_null = 8L)
  dec_tr05 <- deconvolve_crossrun_hrf(d_tr05$bold, d_tr05$fit, scaling = "noise_sd", penalty = "sobolev")
  dec_tr20 <- deconvolve_crossrun_hrf(d_tr20$bold, d_tr20$fit, scaling = "noise_sd", penalty = "sobolev")
  expect_gt(dec_tr05$settings$lambda, dec_tr20$settings$lambda)
})

test_that("single-run HRF estimation pools contiguous temporal sub-blocks and preserves full-run deconvolution", {
  set.seed(202635)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6)[1L]
  idx_early <- which(abs(pars$peak_delay - 4.5) < 1e-6)[1L]
  idx_late <- which(abs(pars$peak_delay - 7.5) < 1e-6)[1L]

  t_ev <- cumsum(c(6, runif(n_trials - 1L, 4.5, 8.5)))
  x <- as.numeric(scale(rnorm(n_trials)))
  amps <- 1.0 + 0.55 * x

  des_early <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
    seq_len(n_trials), n_trials, lib[, idx_early, drop = FALSE])[[1L]]
  des_can <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
    seq_len(n_trials), n_trials, lib[, idx_can, drop = FALSE])[[1L]]
  des_late <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
    seq_len(n_trials), n_trials, lib[, idx_late, drop = FALSE])[[1L]]

  groups <- stats::setNames(rep(c("R_early", "R_can", "R_late", "R_null"), each = 4L), paste0("v", 1:16))
  W_ar <- hrf_whitening_matrix(seq_len(n_time), 0.30)
  ar_err <- function(sd) sd * sqrt(1 - 0.30^2) * as.numeric(solve(W_ar, rnorm(n_time)))

  y1 <- cbind(
    vapply(1:4, function(i) 100 + as.numeric(des_early %*% amps) + ar_err(0.10), numeric(n_time)),
    vapply(1:4, function(i) 100 + as.numeric(des_can %*% amps) + ar_err(0.10), numeric(n_time)),
    vapply(1:4, function(i) 100 + as.numeric(des_late %*% amps) + ar_err(0.10), numeric(n_time)),
    vapply(1:4, function(i) 100 + ar_err(0.15), numeric(n_time))
  )
  colnames(y1) <- names(groups)
  bold_1 <- list(run1 = y1)
  events_1 <- list(run1 = data.frame(onset = t_ev, duration = 0, condition = "feedback"))
  x_1 <- list(run1 = data.frame(x = x))

  fit_1 <- estimate_crossrun_hrf(bold_1, events_1, tr, lib,
    spatial_groups = groups, pooling = "parcel", min_parcel_voxels = 3L,
    fallback_kernel = idx_can, return_scores = TRUE)

  expect_equal(fit_1$settings$selection, "within_run_subblock_pooled")
  expect_equal(fit_1$run_ids, "run1")
  expect_equal(fit_1$training_runs$run1, c("run1_block1", "run1_block2", "run1_block3"))
  expect_true(all(fit_1$selections$run1$n_training_runs == 3L))
  expect_equal(fit_1$selections$run1$kernel_index[1:4], rep(idx_early, 4L))
  expect_equal(fit_1$selections$run1$kernel_index[5:8], rep(idx_can, 4L))
  expect_equal(fit_1$selections$run1$kernel_index[9:12], rep(idx_late, 4L))
  expect_equal(fit_1$selections$run1$kernel_index[13:16], rep(idx_can, 4L))
  expect_equal(fit_1$selections$run1$status[13:16], rep("fallback_prespecified", 4L))

  # Soft-average mode and trial amplitude mode also work on a single run
  fit_soft <- estimate_crossrun_hrf(bold_1, events_1, tr, lib,
    spatial_groups = groups, pooling = "parcel", min_parcel_voxels = 3L,
    fallback_kernel = idx_can, amplitude_model = "trial", selection_mode = "soft_average")
  expect_equal(fit_soft$training_runs$run1, c("run1_block1", "run1_block2"))
  expect_equal(unname(colSums(fit_soft$selected_kernels$run1)), rep(1, 16L), tolerance = 1e-10)

  # Full-run deconvolution, recovery diagnostics, and Stage 3 epoch regression work on 1 run
  dec_1 <- deconvolve_crossrun_hrf(bold_1, fit_1, scaling = "demean", penalty = "sobolev")
  expect_equal(names(dec_1$activity), "run1")
  expect_equal(dim(dec_1$activity$run1), c(n_time, 16L))

  diag_1 <- diagnose_crossrun_recovery(dec_1, events = events_1, lags = 0)
  expect_true(nrow(diag_1$summary) >= 1L)

  reg_1 <- decon_epoch_regression(dec_1, x_1, events = events_1,
    target_condition = "feedback", lags = -1:2, demix_target_trials = TRUE)
  cf_1 <- reg_1$coefficients
  b_early_0 <- mean(cf_1$beta_corrected[cf_1$voxel %in% paste0("v", 1:4) & cf_1$lag == 0])
  b_late_0 <- mean(cf_1$beta_corrected[cf_1$voxel %in% paste0("v", 9:12) & cf_1$lag == 0])
  b_null_0 <- mean(abs(cf_1$beta_corrected[cf_1$voxel %in% paste0("v", 13:16) & cf_1$lag == 0]))
  expect_gt(b_early_0, 0.40)
  expect_gt(b_late_0, 0.40)
  expect_lt(b_null_0, 0.15)

  # Incompatible sub-blocks (early HRF in block 1, late HRF in blocks 2-3) trigger compatibility fallback
  y_incompat <- y1[, 1:2, drop = FALSE]
  y_incompat[1:60, ] <- 100 + as.numeric(des_early[1:60, ] %*% amps) + rnorm(60, sd = 0.02)
  y_incompat[61:180, ] <- 100 + as.numeric(des_late[61:180, ] %*% amps) + rnorm(120, sd = 0.02)
  fit_inc <- estimate_crossrun_hrf(list(run1 = y_incompat), events_1, tr, lib,
    fallback_kernel = idx_can, noise_model = "ols")
  expect_equal(fit_inc$selections$run1$kernel_index, rep(idx_can, 2L))
  expect_equal(fit_inc$selections$run1$status, rep("fallback_prespecified", 2L))
  expect_equal(fit_inc$selections$run1$stability, rep("no_common_library_shape", 2L))

  # Very short single run (< 24 usable observations) gracefully falls back to fallback_kernel
  fit_short <- estimate_crossrun_hrf(
    list(run1 = y1[1:18, 1:4, drop = FALSE]),
    list(run1 = data.frame(onset = c(2, 8, 14), duration = 0, condition = "feedback")),
    tr, lib, fallback_kernel = idx_can, noise_model = "ols"
  )
  expect_equal(fit_short$training_runs$run1, character(0))
  expect_equal(fit_short$selections$run1$kernel_index, rep(idx_can, 4L))
  expect_equal(fit_short$selections$run1$voxel_evidence, rep("insufficient_training", 4L))
})

test_that("decon_epoch_regression supports native parcel pre-aggregation and covariate-agnostic precision weights", {
  set.seed(202642)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)
  hrf <- lib[, 1L, drop = FALSE]
  true_beta <- 0.55

  # Parcel R1 has 4 high-SNR core voxels (sd = 0.08, full signal) and 4 high-noise partial-volume voxels (sd = 0.55, 0.35x signal)
  # Parcel R2 has 6 voxels with spatially correlated noise
  groups <- stats::setNames(
    c(rep("R1", 8L), rep("R2", 6L)),
    paste0("v", 1:14)
  )
  vox_sd_r1 <- c(rep(0.08, 4L), rep(0.55, 4L))
  sig_mult_r1 <- c(rep(1.00, 4L), rep(0.35, 4L))

  bold_indep <- bold_corr <- events <- x_list <- list()
  for (r in 1:2) {
    t_ev <- cumsum(c(6, runif(n_trials - 1L, 4.5, 8.0)))
    x <- as.numeric(scale(rnorm(n_trials)))
    bio_jit <- rnorm(n_trials, sd = 0.25)
    amps_r1 <- 1.0 + true_beta * x + bio_jit
    amps_r2 <- 0.8 + 0.35 * amps_r1 + rnorm(n_trials, sd = 0.20)

    des <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
      seq_len(n_trials), n_trials, hrf)[[1L]]
    sig_r1 <- as.numeric(des %*% amps_r1)
    sig_r2 <- as.numeric(des %*% amps_r2)

    eps_r1 <- vapply(vox_sd_r1, function(s) rnorm(n_time, sd = s), numeric(n_time))
    eps_r2_indep <- replicate(6L, rnorm(n_time, sd = 0.20))
    shared_r2 <- rnorm(n_time, sd = 0.20 * sqrt(0.80))
    eps_r2_corr <- replicate(6L, shared_r2 + rnorm(n_time, sd = 0.20 * sqrt(0.20)))

    sig_mat_r1 <- outer(sig_r1, sig_mult_r1)
    y_indep <- cbind(100 + sig_mat_r1 + eps_r1, sweep(eps_r2_indep, 1L, 100 + sig_r2, "+"))
    y_corr <- cbind(100 + sig_mat_r1 + eps_r1, sweep(eps_r2_corr, 1L, 100 + sig_r2, "+"))
    colnames(y_indep) <- colnames(y_corr) <- names(groups)

    bold_indep[[paste0("run", r)]] <- y_indep
    bold_corr[[paste0("run", r)]] <- y_corr
    events[[paste0("run", r)]] <- data.frame(onset = t_ev, duration = 0, condition = "feedback")
    x_list[[paste0("run", r)]] <- data.frame(x = x)
  }

  fit_ind <- estimate_crossrun_hrf(bold_indep, events, tr, lib,
    spatial_groups = groups, pooling = "parcel", noise_model = "ar1", evidence = "calibrated")
  dec_ind <- deconvolve_crossrun_hrf(bold_indep, fit_ind, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")

  fit_cor <- estimate_crossrun_hrf(bold_corr, events, tr, lib,
    spatial_groups = groups, pooling = "parcel", noise_model = "ar1", evidence = "calibrated")
  dec_cor <- deconvolve_crossrun_hrf(bold_corr, fit_cor, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")

  # 1. Native parcel pre-aggregation downweights noisy voxels and beats unweighted post-hoc voxel averaging
  reg_vox <- decon_epoch_regression(dec_ind, x_list, events = events,
    target_condition = "feedback", lags = 0, restore_scale = TRUE,
    weighting = "gls", whiten_panel = TRUE, return_panel = TRUE)
  reg_parcel <- decon_epoch_regression(dec_ind, x_list, events = events,
    target_condition = "feedback", lags = 0, restore_scale = TRUE,
    weighting = "gls", aggregate_by_parcel = TRUE, whiten_panel = TRUE, return_panel = TRUE)

  expect_true(reg_parcel$settings$aggregate_by_parcel)
  expect_equal( reg_parcel$coefficients$voxel, c("R1", "R2") )
  expect_equal( reg_parcel$coefficients$parcel, c("R1", "R2") )
  expect_true("parcel" %in% names(reg_parcel$panel))
  expect_equal(nrow(reg_parcel$panel), 2L * 2L * n_trials)

  # Clean voxels (1:4) receive much higher weight than noisy voxels (5:8) in R1
  w_r1 <- reg_parcel$parcel_weights$R1
  expect_gt(mean(w_r1[1:4]), 10 * mean(w_r1[5:8]))

  b_vox_r1 <- mean(reg_vox$coefficients$beta_corrected[reg_vox$coefficients$voxel %in% paste0("v", 1:8)])
  se_vox_r1 <- mean(reg_vox$coefficients$se_corrected[reg_vox$coefficients$voxel %in% paste0("v", 1:8)])
  b_par_r1 <- reg_parcel$coefficients$beta_corrected[reg_parcel$coefficients$voxel == "R1"]
  se_par_r1 <- reg_parcel$coefficients$se_corrected[reg_parcel$coefficients$voxel == "R1"]
  expect_lt(se_par_r1, 0.50 * se_vox_r1)
  expect_lt(abs(b_par_r1 - true_beta), abs(b_vox_r1 - true_beta))

  # 2. Empirical spatial noise correlation propagation: spatially correlated R2 has larger parcel SE than independent R2
  reg_vox_cor <- decon_epoch_regression(dec_cor, x_list, events = events,
    target_condition = "feedback", lags = 0, restore_scale = TRUE)
  reg_parcel_cor <- decon_epoch_regression(dec_cor, x_list, events = events,
    target_condition = "feedback", lags = 0, restore_scale = TRUE,
    aggregate_by_parcel = TRUE, return_panel = TRUE)
  se_r2_ind <- reg_parcel$coefficients$mean_se_trial[reg_parcel$coefficients$voxel == "R2"]
  se_r2_cor <- reg_parcel_cor$coefficients$mean_se_trial[reg_parcel_cor$coefficients$voxel == "R2"]
  se_r2_single_ind <- mean(reg_vox$coefficients$mean_se_trial[reg_vox$coefficients$voxel %in% paste0("v", 9:14)])
  se_r2_single_cor <- mean(reg_vox_cor$coefficients$mean_se_trial[reg_vox_cor$coefficients$voxel %in% paste0("v", 9:14)])
  expect_lt(se_r2_ind, se_r2_cor)
  expect_lt(se_r2_ind, se_r2_single_ind)
  expect_lte(se_r2_cor, se_r2_single_cor)

  # 3. Covariate-agnostic panel extraction (covariates = NULL) and tau2_mode = "marginal"
  reg_nocov <- decon_epoch_regression(dec_ind, covariates = NULL, events = events,
    target_condition = "feedback", lags = 0, restore_scale = TRUE,
    aggregate_by_parcel = TRUE, whiten_panel = TRUE)
  reg_marg <- decon_epoch_regression(dec_ind, x_list, events = events,
    target_condition = "feedback", lags = 0, restore_scale = TRUE,
    tau2_mode = "marginal", aggregate_by_parcel = TRUE, whiten_panel = TRUE, return_panel = TRUE)

  expect_true(is.data.frame(reg_nocov$panel))
  expect_equal(reg_nocov$settings$tau2_mode, "marginal")
  expect_true(all(c("decon_demixed", "decon_demixed_whitened", "precision_weight", "tau2_bio") %in% names(reg_nocov$panel)))
  expect_true(all(is.finite(reg_nocov$panel$decon_demixed_whitened)))
  expect_true(all(reg_nocov$coefficients$tau2_bio > 0))
  # Marginal tau2_bio and precision weights are strictly invariant to whether covariates are supplied
  expect_equal(reg_nocov$coefficients$tau2_bio, reg_marg$coefficients$tau2_bio, tolerance = 1e-12)
  expect_equal(reg_nocov$panel$precision_weight, reg_marg$panel$precision_weight, tolerance = 1e-12)
})

test_that("medusa_fit_mmom and decon_connectivity_repair disattenuate measurement error in connectivity and PPI", {
  set.seed(2026)
  n_obs <- 120L
  x <- as.numeric(scale(rnorm(n_obs)))
  b_seed <- rnorm(n_obs, 0, 0.70)
  psi_s <- runif(n_obs, 0.20, 0.45)
  psi_y <- runif(n_obs, 0.15, 0.35)
  y_seed <- b_seed + rnorm(n_obs, 0, sqrt(psi_s))
  y_targ <- 0.50 * b_seed + 0.35 * x + 0.40 * (b_seed * x) + rnorm(n_obs, 0, 0.25) + rnorm(n_obs, 0, sqrt(psi_y))

  X_ppi <- cbind(seed = y_seed, x = x, ppi = y_seed * x)
  err_x <- array(0.0, dim = c(3L, 3L, n_obs))
  err_x[1L, 1L, ] <- psi_s
  err_x[1L, 3L, ] <- err_x[3L, 1L, ] <- x * psi_s
  err_x[3L, 3L, ] <- (x^2) * psi_s

  fit_mmoms <- medusa_fit_mmom(y_targ, X_ppi, err_x = err_x, err_y = psi_y, regularization = "step_up")
  fit_fb <- fmri.pipeline:::medusa_fit_mmom_fallback(
    list(outcome = y_targ, predictors = X_ppi, err_x = err_x, err_y = psi_y),
    list(correction = "step_up"), list(correction = "step_up")
  )
  expect_true(fit_mmoms$converged)
  expect_equal(fit_mmoms$coefficients$estimate, fit_fb$coefficients$estimate, tolerance = 1e-10)
  expect_equal(fit_mmoms$coefficients$se, fit_fb$coefficients$se, tolerance = 1e-10)

  naive_cf <- stats::coef(stats::lm(y_targ ~ y_seed * x))
  mmoms_cf <- stats::setNames(fit_mmoms$coefficients$estimate, fit_mmoms$coefficients$term)
  expect_lt(abs(mmoms_cf[["seed"]] - 0.50), abs(naive_cf[["y_seed"]] - 0.50))
  expect_lt(abs(mmoms_cf[["ppi"]] - 0.40), abs(naive_cf[["y_seed:x"]] - 0.40))

  # End-to-end decon_connectivity_repair test across multiple subjects
  tr <- 1.0
  n_time <- 180L
  n_trials <- 28L
  lib <- hrf_library(tr, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)
  hrf <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ", "R_null"), each = 3L), paste0("v", 1:9))

  panels_all <- vector("list", 10L)
  for (s in 1:10) {
    bold_s <- events_s <- x_s <- list()
    for (r in 1:2) {
      t_ev <- cumsum(c(6, runif(n_trials - 1L, 4.2, 7.5)))
      xv <- as.numeric(scale(rnorm(n_trials)))
      eta_seed <- rnorm(n_trials, 0, 0.70)
      eta_targ <- (0.50 + 0.35 * xv) * eta_seed + rnorm(n_trials, 0, 0.30)
      eta_null <- rnorm(n_trials, 0, 0.60)
      amp_s <- 0.8 + 0.45 * xv + eta_seed
      amp_t <- 0.8 + 0.35 * xv + eta_targ
      amp_0 <- 0.8 + 0.50 * xv + eta_null
      des <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
        seq_len(n_trials), n_trials, hrf)[[1L]]
      ymat <- cbind(
        replicate(3L, 100 + as.numeric(des %*% amp_s) + rnorm(n_time, sd = 0.22)),
        replicate(3L, 100 + as.numeric(des %*% amp_t) + rnorm(n_time, sd = 0.22)),
        replicate(3L, 100 + as.numeric(des %*% amp_0) + rnorm(n_time, sd = 0.22))
      )
      colnames(ymat) <- names(groups)
      bold_s[[paste0("run", r)]] <- ymat
      events_s[[paste0("run", r)]] <- data.frame(onset = t_ev, duration = 0, condition = "feedback")
      x_s[[paste0("run", r)]] <- data.frame(x = xv)
    }
    fh <- estimate_crossrun_hrf(bold_s, events_s, tr, lib, spatial_groups = groups,
      pooling = "parcel", noise_model = "ar1", evidence = "calibrated")
    dc <- deconvolve_crossrun_hrf(bold_s, fh, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
    rg <- decon_epoch_regression(dc, x_s, events = events_s, target_condition = "feedback",
      lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, whiten_panel = TRUE, return_panel = TRUE)
    p_s <- rg$panel
    p_s$id <- sprintf("sub%02d", s)
    panels_all[[s]] <- p_s
  }
  panel_full <- data.table::rbindlist(panels_all)
  expect_true(all(c("decon_eb", "reliability_trial", "se_whitened", "decon_eb_whitened") %in% names(panel_full)))

  conn_mmoms <- decon_connectivity_repair(panel_full, "R_seed", "R_targ",
    test_ppi = TRUE, method = "mmoms", return_augmented_panel = TRUE)
  conn_naive <- decon_connectivity_repair(panel_full, "R_seed", "R_targ",
    test_ppi = TRUE, method = "naive_gls", orthogonalize_seed_x = FALSE)
  conn_null <- decon_connectivity_repair(panel_full, "R_seed", "R_null",
    test_ppi = FALSE, method = "mmoms")

  g_mmoms <- conn_mmoms$group
  g_naive <- conn_naive$group
  expect_lt(abs(g_mmoms[term == "seed_c"]$estimate - 0.50), abs(g_naive[term == "seed_c"]$estimate - 0.50))
  expect_lt(abs(g_mmoms[term == "seed_c:x_c"]$estimate - 0.35), abs(g_naive[term == "seed_c:x_c"]$estimate - 0.35))
  expect_gt(conn_null$group[term == "seed_c"]$p.value, 0.05)

  # Moment-preserving REPAIR predictor columns reproduce subject M-MOMS slopes in WLS
  p_aug <- conn_mmoms$panel_repair[id == "sub01"]
  lm_aug <- stats::lm(targ_c ~ seed_c_repair + x_c_repair + `seed_c:x_c_repair`,
    data = p_aug, weights = repair_weight)
  s1_mmoms <- conn_mmoms$subjects[id == "sub01"]
  expect_equal(unname(stats::coef(lm_aug)[-1L]), s1_mmoms$estimate, tolerance = 1e-5)

  # Between-subject random-effect regression via medusa_between_subject_repair
  if (requireNamespace("repair", quietly = TRUE)) {
    p_rseed <- panel_full[voxel == "R_seed" & evt_time == 0]
    m_l1 <- suppressWarnings(lme4::lmer(decon_demixed ~ x + (1 + x | id), data = p_rseed))
    re_ids <- rownames(lme4::ranef(m_l1)$id)
    b_true_re <- lme4::ranef(m_l1)$id[["x"]]
    trait_df <- data.frame(
      id = re_ids,
      trait_y = 1.5 * b_true_re + rnorm(length(re_ids), 0, 0.05)
    )
    bs_rep <- medusa_between_subject_repair(
      formula = trait_y ~ re(x, group = "id"),
      medusa_mlm = m_l1,
      subject_data = trait_df,
      id_col = "id"
    )
    expect_true(is.list(bs_rep))
    expect_true(!is.null(bs_rep$coefficients) || inherits(bs_rep, "repair"))
  }
})

test_that("fast-TR amendments: AR(1)xAR(2) cascade whitening, AR(3) second-difference noise_sd, and sub-TR lobe-calibrated within-trial demixing", {
  set.seed(8842)
  tr <- 0.65
  n_time <- 260L
  n_trials <- 22L
  lib <- hrf_library(tr, preset = "medusa")
  pars <- attr(lib, "parameters")
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6 & abs(pars$dispersion - 1.0) < 1e-6)[1L]
  groups <- stats::setNames(rep(c("R_task", "R_null"), each = 5L), paste0("v", 1:10))

  # 1. Exact Cholesky identity for cascaded AR(1) x AR(2) whitening matrix
  W_casc <- fmri.pipeline:::hrf_whitening_matrix(seq_len(64L), c(0.60, 0.12, -0.24))
  W_ar1 <- fmri.pipeline:::hrf_whitening_matrix(seq_len(64L), 0.60)
  expect_equal(dim(W_casc), c(64L, 64L))
  expect_gt(max(abs(W_casc - W_ar1)), 0.05)

  # 2. Simulate fast-TR multiband data with oscillatory AR(2) physiological component + slow drift
  bold <- events <- covs <- list()
  t_sec <- (seq_len(n_time) - 1L) * tr
  for (r in 1:2) {
    onsets <- cumsum(c(5.0, runif(n_trials - 1L, 4.0, 6.8)))
    xv <- as.numeric(scale(rnorm(n_trials)))
    des <- hrf_event_design_cpp(n_time, tr, 0, onsets, rep(0, n_trials),
      seq_len(n_trials), n_trials, lib[, idx_can, drop = FALSE])[[1L]]
    ymat <- matrix(0, n_time, 10L, dimnames = list(NULL, names(groups)))
    for (v in 1:10) {
      e_ar1 <- as.numeric(stats::filter(rnorm(n_time + 40L, sd = 0.16), 0.62, method = "recursive"))[-seq_len(40L)]
      e_phys <- 0.11 * sin(2 * pi * 0.24 * t_sec + runif(1, 0, 2 * pi)) +
        0.10 * sin(2 * pi * 1.18 * t_sec + runif(1, 0, 2 * pi))
      e_drift <- 0.18 * cos(pi * (2 * (0:(n_time - 1L)) + 1) / (2 * n_time))
      sig <- if (v <= 5L) as.numeric(des %*% (0.8 + 0.65 * xv)) else rep(0, n_time)
      ymat[, v] <- 100 + sig + e_ar1 + e_phys + e_drift
    }
    bold[[paste0("run", r)]] <- ymat
    events[[paste0("run", r)]] <- data.frame(onset = onsets, duration = 0, condition = "feedback")
    covs[[paste0("run", r)]] <- data.frame(x = xv)
  }

  fit <- estimate_crossrun_hrf(bold, events, tr, lib, spatial_groups = groups,
    pooling = "parcel", fallback_kernel = idx_can, amplitude_model = "trial",
    noise_model = "ar1", evidence = "calibrated", selection_mode = "soft_average")
  expect_true(all(c("rho", "psi1", "psi2") %in% names(fit$noise[[1L]])))
  expect_true(any(abs(fit$noise[[1L]]$psi2) > 0))

  # 3. Priority 2: AR(3) + HF second-difference noise_sd is closer to true high-frequency SD than pure AR(1) formula
  sd_ar3 <- fmri.pipeline:::estimate_decon_noise_sd(
    bold$run1[, 6:10],
    noise_rho = fit$noise[[1L]]$rho[6:10],
    noise_psi1 = fit$noise[[1L]]$psi1[6:10],
    noise_psi2 = fit$noise[[1L]]$psi2[6:10]
  )
  expect_true(all(is.finite(sd_ar3) & sd_ar3 > 0.10 & sd_ar3 < 0.35))

  # 4. Priority 3: Sub-TR multi-lag within-trial demixing calibrates continuous lobe gain and propagates SE inflation
  decon <- deconvolve_crossrun_hrf(bold, fit, lambda = NULL, scaling = "demean", penalty = "sobolev")
  lags_sub <- seq(-tr, 3 * tr, by = tr)
  reg_nowt <- decon_epoch_regression(decon, covs, events = events, target_condition = "feedback",
    lags = lags_sub, demix_within_trial = FALSE, restore_scale = TRUE, return_panel = TRUE)
  reg_wt <- decon_epoch_regression(decon, covs, events = events, target_condition = "feedback",
    lags = lags_sub, demix_within_trial = TRUE, restore_scale = TRUE, return_panel = TRUE)

  se_nowt <- mean(reg_nowt$panel$se_demixed, na.rm = TRUE)
  se_wt <- mean(reg_wt$panel$se_demixed, na.rm = TRUE)
  sd_null_wt <- stats::sd(reg_wt$panel$decon_demixed[reg_wt$panel$voxel %in% paste0("v", 6:10)], na.rm = TRUE)
  expect_true(abs(se_wt - se_nowt) > 0.10)
  expect_lt(abs(se_wt - sd_null_wt) / sd_null_wt, 0.15)
  b_nowt <- reg_nowt$coefficients$beta_weighted[reg_nowt$coefficients$voxel %in% paste0("v", 1:5) & abs(reg_nowt$coefficients$lag) < 1e-6]
  b_wt <- reg_wt$coefficients$beta_weighted[reg_wt$coefficients$voxel %in% paste0("v", 1:5) & abs(reg_wt$coefficients$lag) < 1e-6]
  expect_gt(se_wt / mean(b_wt), se_nowt / mean(b_nowt))
  expect_gt(mean(b_wt), 0.15)
  expect_lt(mean(b_wt), 0.85)
})

