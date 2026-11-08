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
