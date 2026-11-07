test_that("invalid signals cannot change valid estimates or tuning in any backend", {
  set.seed(2011)
  good <- matrix(rnorm(60), ncol = 2L)
  input <- cbind(good[, 1L], 100, good[, 2L], NA_real_, NaN, Inf)
  # Partially missing columns are excluded too, without shortening the clock.
  input[1:29, 4L] <- good[1:29, 1L]

  for (backend in c("cpp", "r", "wiener")) {
    for (tune_by in c("global", "signal")) {
      settings <- list(kernel = c(0, 0.25, 0.5, 0.25),
        lambda_grid = c(0.01, 0.1), backend = backend, tune_by = tune_by,
        prewhiten_gcv = TRUE, return_diagnostics = TRUE)
      reference <- do.call(deconvolve_reglin, c(list(BOLDobs = good), settings))
      expect_warning(
        actual <- do.call(deconvolve_reglin, c(list(BOLDobs = input), settings)),
        "4 constant or nonfinite signals"
      )
      expect_equal(dim(actual$activity), c(30L, 6L))
      expect_equal(actual$activity[, c(1L, 3L)], reference$activity, ignore_attr = TRUE)
      expect_true(all(is.na(actual$activity[, c(2L, 4L, 5L, 6L)])))
      expect_equal(actual$signal_status,
        c("ok", "constant", "ok", "nonfinite", "nonfinite", "nonfinite"))
      expect_equal(actual$tuning, reference$tuning)
      expect_equal(actual$ar1_rho, reference$ar1_rho)
      selected <- if (tune_by == "global") 1L else c(1L, 3L)
      expect_equal(actual$lambda[selected], as.numeric(reference$lambda))
      expect_equal(actual$kernel_index[selected], as.integer(reference$kernel_index))
      if (tune_by == "signal") {
        expect_true(all(is.na(actual$lambda[c(2L, 4L, 5L, 6L)])))
        expect_true(all(is.na(actual$kernel_index[c(2L, 4L, 5L, 6L)])))
      }
      expect_equal(actual$scaled_center[c(1L, 3L)], as.numeric(reference$scaled_center))
      expect_equal(actual$scaled_scale[c(1L, 3L)], as.numeric(reference$scaled_scale))
    }
  }
})

test_that("all invalid signals return missing estimates without fitting", {
  for (backend in c("cpp", "r", "wiener")) {
    for (tune_by in c("global", "signal")) {
      expect_warning(
        result <- deconvolve_reglin(cbind(rep(100, 10L), rep(NA_real_, 10L)),
          kernel = c(0, 0.5, 0.5), backend = backend, tune_by = tune_by,
          normalize = FALSE, trim_kernel = FALSE, prewhiten_gcv = TRUE,
          return_diagnostics = TRUE),
        "2 constant or nonfinite signals"
      )
      expect_equal(dim(result$activity), c(12L, 2L))
      expect_true(all(is.na(result$activity)))
      expect_equal(result$signal_status, c("constant", "nonfinite"))
      expect_equal(nrow(result$tuning), 0L)
      expect_true(all(is.na(result$lambda)))
      expect_equal(length(result$lambda), if (tune_by == "global") 1L else 2L)
      expect_true(is.na(result$ar1_rho))
    }
  }
})

test_that("matrix-only output restores invalid columns with fixed lambda", {
  y <- sin(seq_len(30))
  reference <- deconvolve_reglin(y, kernel = c(0, 0.5, 0.5), lambda = 0.1,
    normalize = FALSE, demean = FALSE)
  expect_warning(
    result <- deconvolve_reglin(cbind(100, y), kernel = c(0, 0.5, 0.5),
      lambda = 0.1, normalize = FALSE, demean = FALSE),
    "1 constant or nonfinite signals"
  )
  expect_true(all(is.na(result[, 1L])))
  expect_equal(result[, 2L], reference[, 1L], ignore_attr = TRUE)
})

test_that("deconvolve_reglin recovers a continuous latent signal", {
  set.seed(1001)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  n_time <- 80L
  n_latent <- n_time + length(kernel) - 1L
  latent <- as.numeric(stats::arima.sim(model = list(ar = 0.6), n = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent

  estimated <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 1e-4,
    penalty = "diff1",
    normalize = FALSE,
    demean = FALSE
  )

  expect_equal(dim(estimated), c(n_time, 1L))
  expect_gt(stats::cor(estimated[, 1L], latent[length(kernel):n_latent]), 0.95)
})

test_that("deconvolve_reglin tunes shifted HRF candidates", {
  set.seed(1002)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  true_lag <- 1L
  true_kernel <- fmri.pipeline:::shift_decon_kernel(true_lag, kernel)
  n_time <- 70L
  n_latent <- n_time + length(kernel) - 1L
  latent <- sin(seq(0, 6 * pi, length.out = n_latent)) + stats::rnorm(n_latent, sd = 0.05)
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, true_kernel) %*% latent

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 1e-4,
    penalty = "diff2",
    hrf_lags = c(-1L, 0L, 1L),
    tune_by = "global",
    normalize = FALSE,
    demean = FALSE,
    return_diagnostics = TRUE
  )

  expect_equal(result$kernel_index, 3L)
  expect_gt(stats::cor(result$activity[, 1L], latent[length(kernel):n_latent]), 0.9)
})

test_that("deconvolve_reglin C++ backend matches R backend", {
  set.seed(1003)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  bold <- matrix(stats::rnorm(90), ncol = 3L)

  cpp_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 0.01,
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    return_diagnostics = TRUE,
    backend = "cpp"
  )

  r_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 0.01,
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    return_diagnostics = TRUE,
    backend = "r"
  )

  expect_equal(cpp_result$activity, r_result$activity, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(cpp_result$lambda, r_result$lambda, tolerance = 1e-12)
  expect_equal(cpp_result$kernel_index, r_result$kernel_index)
})

test_that("deconvolve_reglin prewhitened C++ backend matches R backend", {
  set.seed(1004)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  bold <- matrix(stats::filter(stats::rnorm(120), filter = 0.6, method = "recursive"), ncol = 3L)

  cpp_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = c(0.001, 0.01, 0.1),
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    prewhiten_gcv = TRUE,
    return_diagnostics = TRUE,
    backend = "cpp"
  )

  r_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = c(0.001, 0.01, 0.1),
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    prewhiten_gcv = TRUE,
    return_diagnostics = TRUE,
    backend = "r"
  )

  expect_equal(cpp_result$activity, r_result$activity, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(cpp_result$ar1_rho, r_result$ar1_rho, tolerance = 1e-8)
  expect_true(cpp_result$prewhiten_gcv)
})

test_that("deconvolve_reglin Wiener backend recovers frequency-domain latent signal", {
  set.seed(1005)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 90L
  n_latent <- n_time + length(kernel) - 1L
  latent <- as.numeric(stats::arima.sim(model = list(ar = 0.5), n = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent
  bold <- bold + stats::rnorm(n_time, sd = stats::sd(bold) / 20)

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 1e-4,
    normalize = FALSE,
    demean = FALSE,
    backend = "wiener",
    return_diagnostics = TRUE
  )

  expect_equal(dim(result$activity), c(n_time, 1L))
  expect_equal(result$penalty, "wiener")
  expect_gt(stats::cor(result$activity[, 1L], latent[length(kernel):n_latent]), 0.85)
})

test_that("deconvolve_reglin Wiener backend tunes lambda", {
  set.seed(1006)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 80L
  n_latent <- n_time + length(kernel) - 1L
  latent <- sin(seq(0, 4 * pi, length.out = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent
  bold <- bold + stats::rnorm(n_time, sd = stats::sd(bold) / 4)

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = c(1e-5, 1e-3, 1e-1),
    normalize = FALSE,
    demean = FALSE,
    backend = "wiener",
    return_diagnostics = TRUE
  )

  expect_true(result$lambda %in% c(1e-5, 1e-3, 1e-1))
  expect_equal(nrow(result$tuning), 3L)
})

test_that("deconvolve_reglin PSD Wiener backend recovers latent signal", {
  set.seed(1007)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 90L
  n_latent <- n_time + length(kernel) - 1L
  latent <- as.numeric(stats::arima.sim(model = list(ar = 0.6), n = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent
  bold <- bold + as.numeric(stats::arima.sim(model = list(ar = 0.5), n = n_time)) * stats::sd(bold) / 8

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 1e-2,
    normalize = FALSE,
    demean = FALSE,
    backend = "wiener",
    wiener_shrinkage = "psd",
    return_diagnostics = TRUE
  )

  expect_equal(dim(result$activity), c(n_time, 1L))
  expect_equal(result$wiener_shrinkage, "psd")
  expect_equal(nrow(result$tuning), 1L)
  expect_true(is.na(result$tuning$gcv))
  expect_equal(result$df_method, "conditional_fixed_psd")
  expect_gt(stats::cor(result$activity[, 1L], latent[length(kernel):n_latent]), 0.75)
})

test_that("deconvolve_reglin one-standard-error GCV rule selects smoother C++ candidate", {
  set.seed(1008)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 70L
  bold <- matrix(stats::rnorm(n_time * 12L), nrow = n_time)
  lambda_grid <- c(1e-4, 1e-2, 1)

  min_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "cpp",
    return_diagnostics = TRUE,
    gcv_rule = "min"
  )

  one_se_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "cpp",
    return_diagnostics = TRUE,
    gcv_rule = "1se"
  )

  expect_equal(dim(one_se_result$activity), c(n_time, ncol(bold)))
  expect_true("gcv_se" %in% names(one_se_result$tuning))
  expect_gte(one_se_result$lambda, min_result$lambda)
  expect_equal(one_se_result$gcv_rule, "1se")
})

test_that("deconvolve_reglin one-standard-error GCV rule works for Wiener backend", {
  set.seed(1009)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 70L
  bold <- matrix(stats::rnorm(n_time * 10L), nrow = n_time)
  lambda_grid <- c(1e-4, 1e-2, 1)

  min_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "wiener",
    return_diagnostics = TRUE,
    gcv_rule = "min"
  )

  one_se_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "wiener",
    return_diagnostics = TRUE,
    gcv_rule = "1se"
  )

  expect_equal(dim(one_se_result$activity), c(n_time, ncol(bold)))
  expect_true("gcv_se" %in% names(one_se_result$tuning))
  expect_gte(one_se_result$lambda, min_result$lambda)
  expect_equal(one_se_result$gcv_rule, "1se")
})


test_that("explicit linear scaling has correct units and records the applied transform", {
  y <- cbind(100 + sin(seq_len(40)), 200 + 3 * cos(seq_len(40)))
  transforms <- list(none = y, demean = sweep(y, 2, colMeans(y)),
    zscore = scale(y), psc = sweep(sweep(y, 2, colMeans(y)), 2, colMeans(y) / 100, "/"))
  for (backend in c("cpp", "r", "wiener")) {
    for (mode in names(transforms)) {
      fit <- deconvolve_reglin(y, c(0.1, 0.6, 0.3), lambda = 0.1,
        scaling = mode, backend = backend, return_diagnostics = TRUE)
      ref <- deconvolve_reglin(transforms[[mode]], c(0.1, 0.6, 0.3),
        lambda = 0.1, normalize = FALSE, demean = FALSE, backend = backend)
      expect_equal(fit$activity, ref, ignore_attr = TRUE, tolerance = 1e-9)
      expect_equal(fit$scaling, mode)
      expect_equal(sweep(sweep(y, 2, fit$scaled_center), 2, fit$scaled_scale, "/"),
        transforms[[mode]], ignore_attr = TRUE)
    }
    expect_error(deconvolve_reglin(y - 150, c(0.1, 0.6, 0.3), scaling = "psc", backend = backend),
      "positive baseline")
  }
})

test_that("linear estimates and tuning behave predictably across arbitrary input scales", {
  set.seed(20261106)
  y <- matrix(rnorm(120), 40, 3)
  for (backend in c("cpp", "r", "wiener")) {
    ref <- deconvolve_reglin(y, c(0.1, 0.6, 0.3), lambda_grid = c(0.01, 0.1, 1),
      scaling = "none", backend = backend, return_diagnostics = TRUE)
    for (gain in c(1e-4, 2, 1e6)) {
      scaled <- deconvolve_reglin(y * gain, c(0.1, 0.6, 0.3),
        lambda_grid = c(0.01, 0.1, 1), scaling = "none", backend = backend,
        return_diagnostics = TRUE)
      expect_equal(scaled$activity / gain, ref$activity, ignore_attr = TRUE, tolerance = 1e-8)
      expect_equal(scaled$lambda, ref$lambda)
      standardized <- deconvolve_reglin(y * gain + 100 * gain, c(0.1, 0.6, 0.3),
        lambda = 0.1, backend = backend)
      expect_equal(standardized, deconvolve_reglin(y, c(0.1, 0.6, 0.3),
        lambda = 0.1, backend = backend), ignore_attr = TRUE, tolerance = 1e-8)
    }
  }
})

test_that("Wiener effective df equals the trace of the finite observed smoother", {
  # The dense oracle builds each response to an observed unit impulse, then
  # reconvolves it. It does not use the FFT trace formula under test.
  for (n in c(3L, 8L, 17L, 80L)) {
    for (rho in c(0, 0.6)) {
      kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
      fit_kernel <- fmri.pipeline:::ar1_filter_kernel(kernel, rho)
      m <- n + length(kernel) - 1L
      for (lam in c(0, 0.001, 0.1, 10)) {
        x <- fmri.pipeline:::wiener_deconvolve_matrix(diag(n), fit_kernel, lam, m)
        smoother <- fmri.pipeline:::wiener_reconvolve_matrix(x, fit_kernel, n)
        plan <- fmri.pipeline:::wiener_fft_plan(n, m, fit_kernel)
        filter <- fmri.pipeline:::wiener_inverse_filter(plan$h_fft, lam)
        df <- fmri.pipeline:::wiener_degrees_freedom(filter, fit_kernel, n)
        expect_equal(df, sum(diag(smoother)), tolerance = 1e-11)
      }
    }
  }
  fit <- deconvolve_reglin(sin(1:80), kernel, lambda = 0.1, scaling = "none",
    backend = "wiener", return_diagnostics = TRUE)
  expect_equal(fit$tuning$df, 40.32813, tolerance = 1e-6)
})

test_that("Wiener tuning agrees with independent dense smoother scoring", {
  set.seed(1908)
  y <- matrix(rnorm(36), 12, 3)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  lambdas <- c(0.001, 0.01, 0.1, 1)
  for (rho in c(0, 0.5)) {
    fit_kernel <- fmri.pipeline:::ar1_filter_kernel(kernel, rho)
    y_fit <- fmri.pipeline:::apply_ar1_whitening(y, rho)
    scores <- vapply(lambdas, function(lam) {
      x <- fmri.pipeline:::wiener_deconvolve_matrix(diag(12), fit_kernel, lam, 16)
      smoother <- fmri.pipeline:::wiener_reconvolve_matrix(x, fit_kernel, 12)
      colSums((y_fit - smoother %*% y_fit)^2) / (12 - sum(diag(smoother)))^2
    }, numeric(3))
    for (tune in c("global", "signal")) {
      fit <- deconvolve_reglin(y, kernel, lambda_grid = lambdas, scaling = "none",
        backend = "wiener", prewhiten_gcv = rho != 0, prewhiten_rho = rho,
        tune_by = tune, return_diagnostics = TRUE)
      expect_equal(fit$tuning$gcv, colMeans(scores), tolerance = 1e-11)
      expected <- if (tune == "global") lambdas[which.min(colMeans(scores))] else
        lambdas[apply(scores, 1, which.min)]
      expect_equal(fit$lambda, expected)
    }
  }
})

test_that("adaptive PSD Wiener does not use a fixed-filter trace for automatic tuning", {
  args <- list(BOLDobs = sin(1:30), kernel = c(0.1, 0.6, 0.3), backend = "wiener",
    wiener_shrinkage = "psd")
  expect_error(do.call(deconvolve_reglin, args), "Automatic GCV tuning is not supported")
  expect_error(do.call(deconvolve_reglin, c(args, list(lambda = 0.1, hrf_lags = 0:1))),
    "Automatic GCV tuning is not supported")
  for (tune in c("global", "signal")) {
    fit <- do.call(deconvolve_reglin, c(args, list(lambda = 0.1, tune_by = tune,
      return_diagnostics = TRUE)))
    expect_true(all(is.finite(fit$activity)))
    expect_equal(fit$lambda, 0.1)
    expect_true(is.na(fit$tuning$gcv))
  }
})

test_that("optimized C++ fitting preserves reference activity, df and candidate selection", {
  set.seed(267)
  y <- matrix(rnorm(120), 40, 3)
  for (penalty in c("ridge", "diff1", "diff2")) {
    for (tune in c("global", "signal")) {
      for (rho in c(0, 0.4)) {
        args <- list(BOLDobs = y, kernel = c(0.1, 0.6, 0.3),
          lambda_grid = c(0.001, 0.1, 10), hrf_lags = 0:1, penalty = penalty,
          tune_by = tune, prewhiten_gcv = rho != 0, prewhiten_rho = rho,
          return_diagnostics = TRUE)
        actual <- do.call(deconvolve_reglin, c(args, list(backend = "cpp")))
        oracle <- do.call(deconvolve_reglin, c(args, list(backend = "r")))
        expect_equal(actual$activity, oracle$activity, ignore_attr = TRUE, tolerance = 1e-8)
        expect_equal(actual$tuning, oracle$tuning, tolerance = 1e-8)
        expect_equal(actual$lambda, oracle$lambda, ignore_attr = TRUE)
        expect_equal(actual$kernel_index, oracle$kernel_index, ignore_attr = TRUE)
      }
    }
  }
})

test_that("fixed-parameter fast paths agree with diagnostic fits", {
  y <- cbind(sin(1:30), cos(1:30))
  for (backend in c("cpp", "wiener")) {
    for (tune in c("global", "signal")) {
      for (rule in c("min", "1se")) {
        args <- list(BOLDobs = y, kernel = c(0.1, 0.6, 0.3), lambda = 0.1,
          backend = backend, tune_by = tune, gcv_rule = rule,
          prewhiten_gcv = TRUE, prewhiten_rho = 0.4, trim_kernel = FALSE)
        fast <- do.call(deconvolve_reglin, args)
        full <- do.call(deconvolve_reglin, c(args, list(return_diagnostics = TRUE)))
        expect_equal(fast, full$activity, tolerance = 1e-8)
        expect_equal(as.numeric(full$lambda), rep(0.1, if (tune == "global") 1L else 2L))
      }
    }
  }
})

test_that("banded Cholesky and Takahashi trace match dense R reference across kernels and edge cases", {
  set.seed(4102)
  y <- matrix(rnorm(180), 60, 3)
  kernels_list <- list(
    short = c(0.1, 0.6, 0.3),
    canonical = as.numeric(fmri.pipeline:::spm_hrf(1)$hrf[1:12])
  )
  for (k in kernels_list) {
    for (penalty in c("ridge", "diff1", "diff2")) {
      for (rho in c(0, 0.65)) {
        for (rule in c("min", "1se")) {
          args <- list(BOLDobs = y, kernel = k, lambda_grid = c(1e-4, 1e-2, 1, 10),
            hrf_lags = c(-1L, 0L, 1L), penalty = penalty, tune_by = "global",
            gcv_rule = rule, prewhiten_gcv = rho != 0, prewhiten_rho = rho,
            trim_kernel = FALSE, return_diagnostics = TRUE)
          cpp_fit <- do.call(deconvolve_reglin, c(args, list(backend = "cpp")))
          r_fit <- do.call(deconvolve_reglin, c(args, list(backend = "r")))
          expect_equal(cpp_fit$tuning$df, r_fit$tuning$df, tolerance = 1e-10)
          expect_equal(cpp_fit$tuning$rss, r_fit$tuning$rss, tolerance = 1e-9)
          expect_equal(cpp_fit$tuning$gcv, r_fit$tuning$gcv, tolerance = 1e-9)
          expect_equal(cpp_fit$activity, r_fit$activity, ignore_attr = TRUE, tolerance = 1e-8)
          expect_equal(cpp_fit$lambda, r_fit$lambda)
          expect_equal(cpp_fit$kernel_index, r_fit$kernel_index)
        }
      }
    }
  }
})
