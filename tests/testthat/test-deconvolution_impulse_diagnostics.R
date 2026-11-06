test_that("impulse recovery agrees with an independent dense finite inverse", {
  n <- 24L
  kernel <- c(0.1, 0.6, 0.3)
  events <- data.frame(trial = c("early", "middle", "late"), onset = c(1, 10.5, 22))
  k <- length(kernel)
  h <- fmri.pipeline:::decon_convolution_matrix(n, kernel)
  truth <- matrix(0, n, 3)
  truth[2, 1] <- 1; truth[c(11, 12), 2] <- 0.5; truth[23, 3] <- 1
  bold <- h %*% rbind(matrix(0, k - 1, 3), truth)
  p <- fmri.pipeline:::decon_penalty_crossprod(n + k - 1L, "diff2")
  expected <- solve(crossprod(h) + 0.1 * p + diag(1e-8, n + k - 1L), crossprod(h, bold))[k:(n + k - 1L), ]
  for (backend in c("cpp", "r")) {
    result <- deconvolution_impulse_diagnostics(events, kernel, n, 1, 0.1,
      backend = backend, return_timeseries = TRUE)
    for (j in 1:3) {
      d <- result$time_series[result$time_series$source_trial == events$trial[j], ]
      expect_equal(d$truth, truth[, j])
      expect_equal(d$bold, as.numeric(bold[, j]))
      expect_equal(d$activity, expected[, j], tolerance = 1e-8)
      for (target in 1:3) {
        row <- result$windows[result$windows$source_trial == events$trial[j] &
          result$windows$target_trial == events$trial[target], ]
        values <- approx(0:(n - 1), expected[, j], xout = events$onset[target] + (-3:3))$y
        expect_equal(row$n_valid, sum(is.finite(values)))
        expect_equal(row$window_l1, sum(abs(values), na.rm = TRUE), tolerance = 1e-8)
      }
    }
    expect_equal(nrow(result$windows), 9L)
    expect_equal(result$windows$window_l1_ratio[result$windows$source_trial == result$windows$target_trial], rep(1, 3))
    expect_false(result$settings$attribution)
  }
})

test_that("frozen scaling, centering, whitening and padding reproduce the public inverse", {
  events <- data.frame(trial = "test", onset = 16)
  kernel <- c(0.1, 0.6, 0.3)
  for (backend in c("cpp", "r", "wiener")) {
    result <- deconvolution_impulse_diagnostics(events, kernel, 30, 0.8, 0.1,
      time_offset = 1.6, backend = backend, input_scale = 70, input_demean = TRUE,
      prewhiten_rho = 0.4, tail_padding = 3L, return_timeseries = TRUE)
    y <- result$time_series$bold
    expected <- deconvolve_reglin(c((y - mean(y)) / 70, 0, 0, 0), kernel,
      lambda = 0.1, scaling = "none", backend = backend, prewhiten_gcv = TRUE,
      prewhiten_rho = 0.4)[1:30, 1] * 70
    expect_equal(result$time_series$activity, expected, ignore_attr = TRUE)
    expect_equal(sum(result$time_series$truth), 1)
  }
})

test_that("batching and source selection do not change fixed impulse estimates", {
  events <- data.frame(trial = c("001", "A", "B"), onset = c(5, 9.5, 13))
  for (backend in c("cpp", "wiener")) {
    all <- deconvolution_impulse_diagnostics(events, c(0, 0.5, 0.5), 25, 1, 0.1,
      backend = backend, batch_size = 3, return_timeseries = TRUE)
    subset <- deconvolution_impulse_diagnostics(events, c(0, 0.5, 0.5), 25, 1, 0.1,
      backend = backend, source_trials = c("B", "001"), batch_size = 1,
      return_timeseries = TRUE)
    expect_equal(subset$impulses$source_trial, c("B", "001"))
    for (label in c("B", "001")) {
      expect_equal(subset$time_series$activity[subset$time_series$source_trial == label],
        all$time_series$activity[all$time_series$source_trial == label], tolerance = 1e-10)
    }
    expect_equal(nrow(subset$windows), 6L)
  }
})

test_that("missing, outside and unobserved edge impulses remain in diagnostic manifests", {
  events <- data.frame(trial = c("ok", "missing", "outside", "unobserved"), onset = c(5, NA, -3, 19))
  result <- deconvolution_impulse_diagnostics(events, c(0, 0, 1), 20, 1, 0.1,
    return_timeseries = TRUE)
  expect_equal(result$impulses$status, c("ok", "missing_onset", "outside_run", "unobserved_impulse"))
  expect_true(all(is.na(result$impulses$onset_gain[-1])))
  expect_equal(nrow(result$windows), 16L)
  expect_true(all(result$windows$n_valid[result$windows$source_trial != "ok"] == 0))
  expect_true(all(is.na(result$time_series$activity[result$time_series$source_trial != "ok"])))
})

test_that("generating HRF mismatch and overlapping windows are explicit", {
  events <- data.frame(trial = c("one", "two"), onset = c(8, 8))
  result <- deconvolution_impulse_diagnostics(events, c(0.1, 0.6, 0.3), 25, 1, 0.1,
    forward_kernel = c(0.3, 0.6, 0.1))
  expect_equal(result$settings$forward_kernel, c(0.3, 0.6, 0.1))
  expect_equal(result$windows$window_l1_ratio, rep(1, 4))
  expect_null(result$time_series)
  expect_error(deconvolution_impulse_diagnostics(events, c(0.1, 0.6, 0.3), 25, 1, NULL), "number")
  expect_error(deconvolution_impulse_diagnostics(events, c(0.1, 0.6, 0.3), 25, 1, 0.1,
    source_trials = "unknown"), "subset")
})
