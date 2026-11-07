pr_test_env <- function() {
  env <- new.env(parent = globalenv())
  sys.source(source_tree_file("inst", "validation", "parametric_recovery.R"), env)
  env
}

test_that("parametric generator preserves true trial amplitudes and signed component mass", {
  e <- pr_test_env(); config <- e$pr_scenarios()$dependent_trial_variation
  config$voxels <- 4L; config$runs <- 3L; config$n <- 180L
  data <- e$pr_generate(config, 7611L)
  expect_identical(data, e$pr_generate(config, 7611L))
  for (r in seq_along(data$truth)) {
    m <- nrow(data$trials[[r]])
    expect_equal(data$readout[[r]] %*% data$truth[[r]],
      unname(tail(data$amplitudes[[r]], m)), ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(colSums(data$truth[[r]]), colSums(data$amplitudes[[r]]), ignore_attr = TRUE)
    expect_true(any(data$eta_truth[[r]] != 0))
    expect_true(any(data$noise[[r]] != 0))
    expect_equal(data$bold[[r]], 100 + data$signal[[r]] + data$noise[[r]])
  }
  expect_true(any(data$gain < 0))
  expect_true(any(data$gain > 0))
  # Paired η changes do not alter schedules, covariates, voxel gains or noise.
  config$eta_sd <- 0
  zero <- e$pr_generate(config, 7611L)
  expect_equal(data$trials, zero$trials)
  expect_equal(data$noise, zero$noise)
  expect_equal(data$gain, zero$gain)
  for (r in seq_along(data$signal)) {
    expect_equal(data$signal[[r]] - zero$signal[[r]], data$eta_signal[[r]], tolerance = 1e-12)
    expect_equal(data$truth[[r]] - zero$truth[[r]], data$eta_truth[[r]], tolerance = 1e-12)
  }
})

test_that("task-active nulls retain real activity and population confounds", {
  e <- pr_test_env(); config <- e$pr_scenarios()$confounded_null
  config$eta_sd <- 0
  d <- e$pr_generate(config, 7612L)
  expect_equal(config$beta, 0)
  slopes <- lapply(seq_along(d$trials), function(r) {
    a <- d$readout[[r]] %*% d$truth[[r]]
    expect_true(max(abs(d$signal[[r]])) > .1)
    e$pr_run_slope(a, d$trials[[r]], TRUE)
  })
  expect_equal(unlist(slopes), rep(0, length(unlist(slopes))), tolerance = 1e-12)
  # The unadjusted estimand is confounded even before BOLD/inversion.
  unadjusted <- e$pr_run_slope(d$readout[[1]] %*% d$truth[[1]], d$trials[[1]], FALSE)
  expect_true(any(abs(unadjusted) > .01))
})

test_that("run-slope intervals match an independent calculation and normal null calibration", {
  e <- pr_test_env()
  slopes <- matrix(c(-.3, .1, .2, -.1, .4, 0), ncol = 1)
  got <- e$pr_inference(slopes, 0)
  reference <- stats::t.test(slopes[, 1])
  expect_equal(c(got$lower, got$upper), as.numeric(reference$conf.int))
  expect_equal(got$p, reference$p.value)
  set.seed(7613L)
  null <- e$pr_inference(matrix(rnorm(6L * 10000L), 6L), rep(0, 10000L))
  expect_lt(abs(mean(null$rejected) - .05), .01)
  expect_lt(abs(mean(null$covered) - .95), .01)
  zero <- e$pr_inference(matrix(0, 6L, 2L), c(0, 0))
  expect_equal(zero$p, c(1, 1)); expect_true(all(zero$covered))
})

test_that("truth metrics separate population effects, finite-sample slopes and neural variability", {
  e <- pr_test_env(); d <- e$pr_generate(e$pr_scenarios()$moderate_trial_variation, 7614L)
  recovered <- lapply(seq_along(d$truth), function(r)
    list(activity = d$truth[[r]], eta = d$eta_truth[[r]]))
  got <- e$pr_dataset_metrics(recovered, d)
  expect_equal(unname(got[c("neural_r", "amplitude_r", "eta_r", "eta_gain", "eta_observed_r")]), rep(1, 5))
  expect_equal(unname(got[c("neural_nrmse", "amplitude_rmse", "eta_observed_nrmse",
    "adjusted_sample_slope_rmse")]), rep(0, 4), tolerance = 1e-12)
  # η is not forced orthogonal to RPE: true finite-sample slopes vary.
  expect_gt(got[["adjusted_beta_rmse"]], .001)
  suppressed <- lapply(recovered, function(x) list(activity = .5 * x$activity, eta = .5 * x$eta))
  loss <- e$pr_dataset_metrics(suppressed, d)
  expect_equal(loss[["eta_gain"]], .5)
  expect_gt(loss[["amplitude_rmse"]], 0)
})

test_that("continuous and discrete forward references agree at integer event times", {
  e <- pr_test_env(); n <- 160L; times <- c(20, 60, 110)
  h <- e$pr_hrf(0:32, 5)
  continuous <- vapply(times, function(at) e$pr_hrf(0:(n - 1L) - at, 5), numeric(n))
  expect_equal(e$pr_discrete_forward(n, h) %*% e$pr_deposit(n, times), continuous)
  expect_equal(sum(h / sum(h)), 1)
  fractional <- vapply(times + .5, function(at) e$pr_hrf(0:(n - 1L) - at, 5), numeric(n))
  expect_gt(max(abs(e$pr_discrete_forward(n, h) %*% e$pr_deposit(n, times + .5) - fractional)), .001)
})

test_that("cached grouped simulation inverse matches the public frozen fit", {
  e <- pr_test_env(); config <- e$pr_scenarios()$moderate_trial_variation
  config$voxels <- 4L; config$runs <- 3L; config$n <- 180L
  d <- e$pr_generate(config, 7615L)
  fit <- estimate_crossrun_hrf(d$bold, d$events, 1, hrf_library(1),
    spatial_groups = d$groups, spatial_strength = .02, amplitude_model = "trial")
  public <- deconvolve_crossrun_hrf(d$bold, fit, .01, scaling = "demean", backend = "cpp")
  cache <- new.env(parent = emptyenv())
  for (r in seq_along(d$bold)) {
    h <- fit$kernels[, fit$selections[[r]]$kernel_index, drop = FALSE]
    got <- e$pr_apply(d$bold[[r]], d$eta_signal[[r]], h, fit$noise[[r]]$rho, .01, cache)
    expect_equal(got$activity, public$activity[[r]], ignore_attr = TRUE, tolerance = 1e-8)
    expect_equal(e$pr_apply(d$signal[[r]], d$eta_signal[[r]], h, fit$noise[[r]]$rho, .01, cache)$eta,
      got$eta, tolerance = 1e-12)
  }
  count <- length(ls(cache))
  e$pr_apply(d$bold[[1]], d$eta_signal[[1]], h, fit$noise[[3]]$rho, .01, cache)
  expect_lte(length(ls(cache)), count)
})

test_that("Monte Carlo uncertainty uses datasets as independent units", {
  e <- pr_test_env()
  x <- data.frame(scenario = "null", replicate = 1:4, method = "test", lambda = .01,
    adjusted_false_positive = c(0, .1, .2, .1))
  result <- e$pr_summarize(x)
  expect_equal(result$n, 4)
  expect_equal(result$mean, .1)
  expect_equal(result$mc_se, sd(x$adjusted_false_positive) / 2)
  expect_equal(result$lower, max(0, .1 - qt(.975, 3) * result$mc_se))
  x$adjusted_false_positive[] <- 0
  boundary <- e$pr_summarize(x)
  expect_equal(boundary$lower, 0)
  expect_equal(boundary$upper, 1 - .05^(1 / 4))
})

test_that("paired diagnostic references freeze fitting and change only their named inputs", {
  withr::local_dir(pkg_root)
  e <- new.env(parent = globalenv())
  sys.source(file.path("inst", "validation", "run_parametric_recovery.R"), e)
  config <- e$pr_scenarios()$moderate_trial_variation
  config$voxels <- 4L; config$runs <- 3L; config$n <- 180L
  data <- e$pr_generate(config, 7616L); fit <- e$pr_fit(data)
  cache <- new.env(parent = emptyenv())
  learned <- e$pr_recover(data, fit, "learned_estimated_noise", .01, cache)
  clean <- e$pr_recover(data, fit, "learned_frozen_noiseless", .01, cache)
  fixed <- e$pr_recover(data, fit, "learned_frozen_known_noise", .01, cache)
  for (r in seq_along(learned)) {
    h <- fit$kernels[, fit$selections[[r]]$kernel_index, drop = FALSE]
    injected <- e$pr_apply(data$noise[[r]], data$eta_signal[[r]], h, fit$noise[[r]]$rho, .01, cache)
    expect_equal(learned[[r]]$activity - clean[[r]]$activity, injected$activity, tolerance = 1e-9)
    expect_equal(learned[[r]]$eta, clean[[r]]$eta)
    reference <- e$pr_apply(data$signal[[r]] + data$noise[[r]], data$eta_signal[[r]],
      h, data$config$noise_rho, .01, cache)
    expect_equal(fixed[[r]], reference)
  }
  d <- data.frame(scenario = "test", replicate = rep(1:4, 2),
    method = rep(c("learned_estimated_noise", "oracle_hrf_estimated_noise"), each = 4),
    lambda = .01, amplitude_rmse = c(.1, .2, .3, .4, .05, .1, .15, .2))
  result <- e$pr_paired(d)
  row <- result[result$method == "shape_estimation", ]
  expect_equal(row$n, 4L); expect_equal(row$mean, .125)
  d$adjusted_coverage <- c(rep(.5, 4), rep(1, 4))
  row <- subset(e$pr_paired(d), method == "shape_estimation" & metric == "adjusted_coverage")
  expect_equal(row$lower, -.5); expect_equal(row$upper, -.5)
  # A completed/frozen experiment cannot silently be overwritten.
  tmp <- tempfile(); dir.create(tmp); writeLines("frozen", file.path(tmp, "settings.rds"))
  expect_error(e$run_parametric_recovery(tmp, "development", 1L, 1L), "empty output directory")
})

test_that("report assembly rejects incomplete, duplicated, or changed evaluation provenance", {
  withr::local_dir(pkg_root)
  e <- new.env(parent = globalenv())
  sys.source(file.path("inst", "validation", "report_parametric_recovery.R"), e)
  tmp <- tempfile(); dir.create(tmp)
  scripts <- file.path("inst", "validation", c("parametric_recovery.R", "run_parametric_recovery.R"))
  settings <- list(phase = "evaluation", scenarios = e$pr_scenarios()["active_null"],
    counts = c(active_null = 1L), lambdas = c(.001, .01), methods = c("a", "b"),
    seed_base = 18000000L, script_md5 = tools::md5sum(scripts))
  saveRDS(settings, file.path(tmp, "settings.rds"))
  metrics <- expand.grid(scenario = "active_null", replicate = 1L, method = settings$methods,
    lambda = settings$lambdas, stringsAsFactors = FALSE)
  metrics$amplitude_rmse <- 1
  write.csv(metrics, file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  diagnostic <- data.frame(scenario = "active_null", replicate = 1L, seed = 18050001L)
  write.csv(diagnostic, file.path(tmp, "fit_diagnostics.csv"), row.names = FALSE)
  expect_equal(nrow(e$pr_collect(tmp)$metrics), 4L)
  write.csv(rbind(metrics, metrics[1, ]), file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  expect_error(e$pr_collect(tmp), "Incomplete or duplicate")
  write.csv(metrics[-1, ], file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  expect_error(e$pr_collect(tmp), "Incomplete or duplicate")
  write.csv(metrics, file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  diagnostic$seed <- 18050002L
  write.csv(diagnostic, file.path(tmp, "fit_diagnostics.csv"), row.names = FALSE)
  expect_error(e$pr_collect(tmp), "Incorrect dataset seed")
  settings$script_md5[] <- "changed"
  saveRDS(settings, file.path(tmp, "settings.rds"))
  expect_error(e$pr_collect(tmp), "sources changed")
})
