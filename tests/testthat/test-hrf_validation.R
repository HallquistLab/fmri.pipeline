test_that("independent continuous generator agrees with gamma integrals and sampled shapes", {
  env <- new.env(parent = globalenv())
  sys.source(source_tree_file("inst", "validation", "hrf_recovery.R"), env)
  for (delay in c(5, 6, 7)) {
    h <- env$hrfv_response(0:32, delay, 1)
    expected <- hrf_library(1, peak_delays = delay, dispersions = 1)
    expect_equal(h / sum(h), as.numeric(expected[, 1]), tolerance = 1e-12)
    integral <- integrate(function(t) env$hrfv_response(t, delay, 1), 0, 32)
    expect_equal(integral$value, 1, tolerance = 1e-8)
  }
  for (at in c(-1, 1, 5, 10, 20, 33)) {
    lo <- max(0, at - 32); hi <- min(2.5, at)
    quadrature <- if (hi > lo) integrate(function(s) env$hrfv_response(at - s, 6, 1), lo, hi,
      subdivisions = 200L, rel.tol = 1e-9)$value / 2.5 else 0
    expect_equal(env$hrfv_response(at, 6, 1, duration = 2.5), quadrature, tolerance = 1e-6)
  }
})

test_that("resolution-limited truth preserves mass and does not enter forward generation", {
  env <- new.env(parent = globalenv())
  sys.source(source_tree_file("inst", "validation", "hrf_recovery.R"), env)
  truth <- env$hrfv_project_truth(100, c(10.25, 40.7), c(-2, 3), c(0, 2.5))
  expect_equal(sum(truth), 1, tolerance = 1e-12)
  expect_equal(truth[11:12], c(-1.5, -.5))
  config <- env$hrfv_scenarios()$shifted
  first <- env$hrfv_generate(config, 817)
  env$hrfv_project_truth <- function(n, ...) rep(0, n)
  second <- env$hrfv_generate(config, 817)
  expect_identical(first$bold, second$bold)
  expect_true(all(second$truth[[1]] == 0))
  expect_true(any(first$truth[[1]] != 0))
})

test_that("recovery metrics recognize exact truth without hiding amplitude errors", {
  env <- new.env(parent = globalenv())
  sys.source(source_tree_file("inst", "validation", "hrf_recovery.R"), env)
  data <- env$hrfv_generate(env$hrfv_scenarios()$shifted, 817)
  x <- data$truth[[1]]
  exact <- env$hrfv_metrics(x, x, data$trials[[1]], data$kernels[[1]], data$kernels[[1]])
  expect_equal(unname(exact[c("neural_r", "rpe_gain", "rpe_area_gain")]), rep(1, 3))
  expect_equal(unname(exact[c("neural_nrmse", "hrf_nmse", "epoch_peak_error",
    "neighbor_leakage", "rpe_gain_error", "rpe_area_gain_error")]), rep(0, 6), tolerance = 1e-12)
  attenuated <- env$hrfv_metrics(.5 * x, x, data$trials[[1]], data$kernels[[1]], data$kernels[[1]])
  expect_equal(unname(attenuated[c("neural_r", "rpe_gain", "rpe_area_gain", "neural_nrmse")]),
    c(1, .5, .5, .5), tolerance = 1e-12)
  expect_equal(unname(attenuated[c("rpe_gain_error", "rpe_area_gain_error")]), rep(.5, 2))
})

test_that("Monte Carlo comparisons pair datasets and do not count runs or voxels as replicates", {
  env <- new.env(parent = globalenv())
  sys.source(source_tree_file("inst", "validation", "hrf_recovery.R"), env)
  d <- data.frame(scenario = "test", replicate = rep(1:4, 2),
    method = rep(c("canonical", "learned"), each = 4), lambda = .01,
    neural_r = c(.1, .2, .3, .4, .15, .3, .45, .6))
  result <- env$hrfv_summary(d)
  expect_true(all(result$summary$n == 4L))
  expect_equal(result$paired$n, 4L)
  expect_equal(result$paired$delta, .125)
  expect_equal(result$paired$mc_se, sd(c(.05, .1, .15, .2)) / 2)
  expect_equal(result$paired$improved_fraction, 1)
  expect_equal(result$paired$lower, .125 - qt(.975, 3) * sd(c(.05, .1, .15, .2)) / 2)
})

test_that("cached Monte Carlo application agrees with the public inverse interface", {
  env <- new.env(parent = globalenv())
  sys.source(source_tree_file("inst", "validation", "hrf_recovery.R"), env)
  data <- env$hrfv_generate(env$hrfv_scenarios()$two_runs, 818)
  kernels <- hrf_library(1)
  fit <- estimate_crossrun_hrf(data$bold, data$events, 1, kernels,
    spatial_groups = data$groups, spatial_strength = .02, amplitude_model = "trial", noise_model = "ols", evidence = "legacy")
  public <- deconvolve_crossrun_hrf(data$bold, fit, .01, scaling = "demean")
  cache <- new.env(parent = emptyenv())
  for (r in seq_along(data$bold)) {
    selected <- fit$kernels[, fit$selections[[r]]$kernel_index, drop = FALSE]
    got <- env$hrfv_apply(data$bold[[r]], selected, .01, cache)
    expect_equal(as.vector(got), as.vector(public$activity[[r]]), tolerance = 1e-9)
  }
})
